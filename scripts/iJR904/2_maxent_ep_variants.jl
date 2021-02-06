import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    using Serialization

    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

    import Chemostat
    import Chemostat.LP: MathProgBase
    const Ch = ChK.Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads
end

## -------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_dat_"

const INDEX = UJL.DictTree()
function dat_file(name; kwargs...)
    fname = UJL.mysavename(name, "jls"; kwargs...)
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
end

## -------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED

## -------------------------------------------------------------------
# EXPECTED
let
    # global setup
    method = EXPECTED

    # orig model
    iterator = Kd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator

        # prepare model
        model = base_model(exp)
        objidx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
        M, N = size(model)
        exp_growth = Kd.val(:D, exp)
        growth_ub = ChU.ub(model, iJR.KAYSER_BIOMASS_IDER)
        feasible = exp_growth < growth_ub
        ChU.ub!(model, iJR.KAYSER_BIOMASS_IDER, growth_ub * 1.1) # open a beat helps EP

        lock(WLOCK) do
            nzabs_range = ChU.nzabs_range(model.S)
            @info("Starting... ", 
                exp, method,
                size(model), nzabs_range, 
                feasible,
                threadid()
            ); println()
        end
        # !feasible && continue

        # simulation
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        dat = isfile(datfile) ? deserialize(datfile) : Dict()
        epouts = get!(dat, :epouts, Dict())
        init_len = length(epouts)
        beta_vec = zeros(N)
        approach_status = get!(dat, :approach_status, :running)
        if approach_status == :finished 
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                @show INDEX[method, :DFILE, exp]
            end
            continue
        end
        convth = 0.05
        
        # log approach
        epout_seed = isempty(epouts) ? nothing : epouts[maximum(keys(epouts))]
        betas = [0.0; 10.0.^(3:0.05:15)]
        nan_beta = first(betas)

        for approach in [:log_approach, :linear_approach]
            
            lock(WLOCK) do
                @info("Starting", exp, approach, 
                    length(epouts),
                    threadid()
                ); println()
            end

            for beta in betas

                nan_beta = beta
                haskey(epouts, beta) && continue

                beta_vec[objidx] = beta
                epout = nothing
                try
                    epout = ChEP.maxent_ep(model; 
                        beta_vec, alpha = Inf, damp = 0.9, epsconv = 1e-4, 
                        maxvar = 1e50, minvar = 1e-50, verbose = false, solution = epout_seed,
                        maxiter = 1000
                    )
                catch err; end

                # info
                ep_growth = isnothing(epout) ? 0.0 : ChU.av(model, epout, objidx)
                lock(WLOCK) do
                    @info("Results", exp, beta, 
                        exp_growth, growth_ub, ep_growth, 
                        length(epouts),
                        threadid()
                    ); println()
                end

                # error conditions
                fail = isnothing(epout) || isnan(ep_growth) || ep_growth == 0.0 
                fail && break

                # updating
                epout_seed = epouts[beta] = epout

                # convergence
                converr = abs(ep_growth - exp_growth)/exp_growth
                conv = converr < convth || ep_growth > exp_growth 
                conv && break

            end # for betas

            # Catching
            # update = init_len != length(epouts)
            update = true
            update && lock(WLOCK) do
                serialize(datfile, dat)
                @info("Catching", exp,  
                    length(epouts),
                    basename(datfile),
                    threadid()
                ); println()
            end

            # lineal approach
            last_beta = maximum(keys(epouts))
            step = abs(last_beta - nan_beta) / 100.0
            iszero(step) && break
            betas = range(last_beta, 1.0e15; step)

        end # for approach

        # saving
        lock(WLOCK) do
            INDEX[method, :DFILE, exp] = datfile
            dat[:approach_status] = :finished
            dat[:model] = model |> ChU.compressed_model
            serialize(datfile, dat)
            @info("Finished", exp,  
                length(epouts),
                threadid()
            ); println()
        end
       
    end # for (exp, D)
end


## -------------------------------------------------------------------
# Further convergence
let
    method = EXPECTED
    objider = iJR.KAYSER_BIOMASS_IDER

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        model, epouts = ChU.uncompressed_model(dat[:model]) , dat[:epouts]
        
        exp_growth = Kd.val(:D, exp)
        exp_beta = maximum(keys(epouts))
        exp_epout = epouts[exp_beta]

        lock(WLOCK) do
            @info("Converging...", exp, method,
                exp_beta, exp_epout.status, 
                threadid()
            ); println()
        end
        converg_status = get!(dat, :converg_status, :undone)
        converg_status == :done && continue

        new_epout = nothing
        if exp_epout.status == :unconverged
            try;
                objidx = ChU.rxnindex(model, objider)
                beta_vec = zeros(size(model, 2)); 
                beta_vec[objidx] = exp_beta
                new_epout = ChEP.maxent_ep(model; 
                    beta_vec, alpha = Inf, damp = 0.98, epsconv = 1e-4, 
                    maxvar = 1e50, minvar = 1e-50, verbose = false, 
                    solution = exp_epout, maxiter = 5000
                )
            catch err; @warn("ERROR", err); println() end

            ep_growth = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, objider)
            fail = isnan(ep_growth) || ep_growth == 0.0 
            epouts[exp_beta] = fail ? exp_epout : new_epout
        end
        
        # Saving
        lock(WLOCK) do
            @info("Saving...", exp, exp_beta, 
                exp_epout.status,
                threadid()
            ); println()
        end
        dat[:converg_status] = :done
        serialize(datfile, dat)
    end
end

## -------------------------------------------------------------------
# BOUNDED
let
    method = BOUNDED
    objider = iJR.KAYSER_BIOMASS_IDER
    costider = iJR.COST_IDER
    biomass_f = 0.01

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        if isfile(datfile)
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                @info("Cached loaded (skipping)",
                    exp, D, datfile, threadid()
                ); println()
            end
            continue
        end

        # setup
        model = base_model(exp)
        objidx = ChU.rxnindex(model, objider)
        M, N = size(model)
        exp_growth = Kd.val("D", exp)
        fbaout = ChLP.fba(model, objider, costider)
        fba_growth = ChU.av(model, fbaout, objider)
        ub_growth = min(fba_growth, exp_growth)
        ChU.ub!(model, objider, ub_growth * (1.0 + biomass_f))
        ChU.lb!(model, objider, ub_growth * (1.0 - biomass_f))
        model = ChLP.fva_preprocess(model, 
            check_obj = objider,
            verbose = false
        )

        lock(WLOCK) do
            @info("Doing... ", 
                exp, method, 
                D, threadid()
            ); println()
        end

        # maxent
        epout = ChEP.maxent_ep(model; 
            alpha = Inf, damp = 0.985, epsconv = 1e-4, 
            verbose = false, maxiter = 5000
        )
            
        # storing
        lock(WLOCK) do
            # Storing
            dat = Dict()
            dat[:epouts] = Dict(0.0 => epout)
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)
            INDEX[method, :DFILE, exp] = datfile

            @info("Finished ", exp, threadid())
            println()
        end
    end
end

## -------------------------------------------------------------------
# HOMO
let
    method = HOMO
    objider = iJR.KAYSER_BIOMASS_IDER

    iterator = Kd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        if isfile(datfile)
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                @info("Cached loaded (skipping)",
                    exp, D, datfile, threadid()
                ); println()
            end
            continue
        end

        # setup
        model = base_model(exp)
        objidx = ChU.rxnindex(model, objider)

        lock(WLOCK) do
            @info("Doing... ", 
                exp, method, 
                D, threadid()
            ); println()
        end

        # maxent
        epout = ChEP.maxent_ep(model; 
            alpha = Inf, damp = 0.985, epsconv = 1e-4, 
            verbose = false, maxiter = 5000
        )
            
        # storing
        lock(WLOCK) do
            # Storing
            dat = Dict()
            dat[:epouts] = Dict(0.0 => epout)
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)
            INDEX[method, :DFILE, exp] = datfile

            @info("Finished ", exp, threadid())
            println()
        end
    end
end

## -------------------------------------------------------------------
# save index
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = false)