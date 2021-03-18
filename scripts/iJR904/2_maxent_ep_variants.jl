import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    using Serialization

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    const iJR = ChF.iJR904
    const Fd = ChF.FolsomData # experimental data
    const Bd = ChF.BegData    # cost data

    import Chemostat
    import Chemostat.LP: MathProgBase
    const Ch = ChF.Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads
    UJL.set_cache_dir(iJR.MODEL_CACHE_DIR)
end

## ----------------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_dat_"

const INDEX = UJL.DictTree()
function dat_file(name; kwargs...)
    fname = UJL.mysavename(name, "jls"; kwargs...)
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
end

## ----------------------------------------------------------------------------
function load_model(exp::Int, modelkey::String = "fva_models")
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS[modelkey][exp]
    model |> ChU.uncompressed_model
end

function load_model(modelkey::String = "max_model")
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS[modelkey]
    model |> ChU.uncompressed_model
end

## ----------------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_MAX_POL                = :ME_MAX_POL                 # 
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
function check_cache(datfile, method)
    thid = threadid()
    if isfile(datfile)
        lock(WLOCK) do
            INDEX[method, :DFILE, exp] = datfile
            @info("Cached loaded (skipping)",
                exp, datfile, thid
            )
            println()
        end
        return true
    end
    return false
end

## -------------------------------------------------------------------
# ME_MAX_POL
let
    method = ME_MAX_POL

    # Feed jobs
    Ch = Channel(nthreads()) do ch
        cGLCs = Fd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    # @threads for thid in 1:nthreads()
    for thid in 1:1 # Test
        for (exp, cGLC) in Ch
            
            ## -------------------------------------------------------------------
            # handle cache
            datfile = dat_file(string(DAT_FILE_PREFFIX, method); exp)
            check_cache(datfile, method) && continue

            ## -------------------------------------------------------------------
            # SetUp
            # @info("Set up") # Test
            model =  load_model("max_model")
            # model =  load_model(exp, "fva_models")
            M, N = size(model)
            biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            glcidx = ChU.rxnindex(model, iJR.GLC_EX_IDER)
            exp_growth = Fd.val("D", exp)
            biom_lb, biom_ub = ChU.bounds(model, iJR.BIOMASS_IDER)
            if biom_ub < exp_growth
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = :unfeasible
                    @info("Not feasible (skipping)", 
                        biom_ub ,exp_growth, 
                        thid
                    ); println()
                end
                continue
            end

            cgD_X = -Fd.cval(:GLC, exp) * Fd.val(:D, exp) / Fd.val(:X, exp)
            biom_beta = 0.0
            biom_betas = [biom_beta]
            vg_beta = 0.0
            vg_betas = [vg_beta]
            vg_avPME = 0.0
            vg_avPME_vgb0 = 0.0
            biom_avPME = 0.0
            biom_avPME_vgb0 = 0.0
            biom_diff = 0.0
            vg_diff = 0.0
            beta_vec = zeros(N)
            epouts = Dict()
            epout = nothing
            epout_vgb0 = nothing

            gdth = 1e-3
            upfrec_time = -1 # Test
            stw = 10
            stth = 0.1

            rounditer = 1
            while true
                ## -------------------------------------------------------------------
                # Z GRAD DESCEND: Match biomass momentums
                let
                    target = exp_growth
                    x0 = biom_beta
                    maxΔx = max(abs(biom_beta) * 0.05, 1e4)
                    x1 = x0 + maxΔx * 0.01
                    
                    last_uptime = time()
                    gdit = 1

                    ## -------------------------------------------------------------------
                    function z_fun(gdmodel)

                        biom_beta = UJL.gd_value(gdmodel)
            
                        beta_vec[biomidx] = biom_beta
                        beta_vec[glcidx] = vg_beta
                    
                        # @info("Z MaxEnt") # Test
                        epout = ChEP.maxent_ep(model; 
                            beta_vec,
                            alpha = Inf,
                            maxiter = 5000,  
                            epsconv = 1e-3, 
                            # epsconv = 1e-2, # Test
                            # verbose = false, 
                            verbose = true, # Test
                            solution = epout
                        )

                        biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
                        vg_avPME = ChU.av(model, epout, iJR.GLC_EX_IDER)
                        biom_diff = abs(biom_avPME - exp_growth)
                        vg_diff = abs(vg_avPME - cgD_X)

                        update = gdit == 1 || abs(last_uptime - time()) > upfrec_time || 
                            epout.status != ChEP.CONVERGED_STATUS
            
                        update && lock(WLOCK) do
                            @info(
                                "z grad descent... ", 
                                exp, rounditer, gdit, 
                                epout.status, epout.iter, 
                                (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                                (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                                (biom_beta, vg_beta), 
                                thid
                            ); println()
                            last_uptime = time()
                        end
                        
                        gdit += 1
                        return biom_avPME
                    end

                    ## -------------------------------------------------------------------
                    gdmodel = UJL.grad_desc(z_fun; 
                        x0, x1, gdth, maxΔx, 
                        target, maxiter = 5000, 
                        verbose = false
                    )
                    biom_beta = UJL.gd_value(gdmodel)
                end

                ## -------------------------------------------------------------------
                # CHECK VG VALIDITY
                let
                    beta_vec[biomidx] = biom_beta
                    beta_vec[glcidx] = 0.0

                    # @info("Z b0 MaxEnt") # Test
                    epout_vgb0 = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        maxiter = 5000,  
                        epsconv = 1e-3, 
                        # epsconv = 1e-2, # Test
                        # verbose = false, 
                        verbose = true, # Test
                        solution = epout
                    )   
                    biom_avPME_vgb0 = ChU.av(model, epout_vgb0, biomidx)
                    vg_avPME_vgb0 = ChU.av(model, epout_vgb0, glcidx)
                end

                ## -------------------------------------------------------------------
                # Force vg boundary
                let
                    if abs(vg_avPME_vgb0) <= abs(cgD_X)
                        vg_beta = 0.0
                        epout = epout_vgb0
                    else
                        ## -------------------------------------------------------------------
                        # VG GRAD DESCEND: Match biomass momentums
                        target = cgD_X * 0.95
                        x0 = vg_beta
                        maxΔx = max(abs(vg_beta) * 0.05, 1e2)
                        x1 = x0 + maxΔx * 0.1
                
                        last_uptime = time()
                        gdit = 1
        
                        ## -------------------------------------------------------------------
                        function vg_fun(gdmodel)
        
                            vg_beta = UJL.gd_value(gdmodel)
                
                            beta_vec[biomidx] = biom_beta
                            beta_vec[glcidx] = vg_beta
                        
                            # @info("Vg b0 MaxEnt") # Test
                            epout = ChEP.maxent_ep(model; 
                                beta_vec,
                                alpha = Inf,
                                maxiter = 5000,  
                                epsconv = 1e-3, 
                                # epsconv = 1e-2, # Test
                                # verbose = false, 
                                verbose = true, 
                                solution = epout
                            )
        
                            biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
                            vg_avPME = ChU.av(model, epout, iJR.GLC_EX_IDER)
                            biom_diff = abs(biom_avPME - exp_growth)
                            vg_diff = abs(vg_avPME - cgD_X)

                            update = gdit == 1 || abs(last_uptime - time()) > upfrec_time || 
                                epout.status != ChEP.CONVERGED_STATUS
                
                            update && lock(WLOCK) do
                                @info(
                                    "vg grad descent... ", 
                                    exp, rounditer, gdit, 
                                    epout.status, epout.iter, 
                                    (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                                    (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                                    (biom_beta, vg_beta), 
                                    thid
                                ); println()
                                last_uptime = time()
                            end
                            
                            gdit += 1
                            return vg_avPME
                        end

                        ## -------------------------------------------------------------------
                        gdmodel = UJL.grad_desc(vg_fun; 
                            x0, x1, gdth, maxΔx, 
                            target, maxiter = 5000, 
                            verbose = false
                        )
                        vg_beta = UJL.gd_value(gdmodel)

                    end # if abs(vg_avPME_vgb0) <= abs(cgD_X)
                end

                ## -------------------------------------------------------------------
                # collect betas
                push!(biom_betas, biom_beta)
                push!(vg_betas, vg_beta)
                epouts[(biom_beta, vg_beta)] = epout

                ## -------------------------------------------------------------------
                # convergence
                hasvalid_moments = (abs(vg_avPME_vgb0) <= abs(cgD_X) && 
                    abs(biom_avPME_vgb0 - exp_growth)/exp_growth < gdth)
                isbeta_stationary = UJL.is_stationary(biom_betas, stth, stw) && 
                    UJL.is_stationary(vg_betas, stth, stw)
                conv = hasvalid_moments && isbeta_stationary

                ## -------------------------------------------------------------------
                lock(WLOCK) do
                    @info("Round Done", 
                        exp, rounditer, conv,
                        epout.status, epout.iter, 
                        (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                        (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                        (biom_beta, vg_beta), 
                        thid
                    ); println()
                end
                
                conv && break
                rounditer += 1

                # Test
                break
            end

            ## -------------------------------------------------------------------
            lock(WLOCK) do

                # Storing
                dat = Dict()
                dat[:exp_beta] = (biom_beta, vg_beta)
                dat[:epouts] = epouts
                dat[:model] = model |> ChU.compressed_model

                # caching
                # serialize(datfile, dat) # Test
                INDEX[method, :DFILE, exp] = datfile

                @info("Finished ",
                    exp, rounditer,  
                    epout.status, epout.iter, 
                    (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                    (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                    (biom_beta, vg_beta), 
                    thid
                ); println()
            end

            break # Test

        end # for exp, cGLC
    end # for thid

end

## -------------------------------------------------------------------
# # ME_Z_EXPECTED_G_MOVING
# let
#     method = ME_Z_EXPECTED_G_MOVING

#     # Feed jobs
#     Ch = Channel(1) do ch
#         cGLCs = Fd.val("cGLC")
#         for (exp, cGLC)  in enumerate(cGLCs)
#             put!(ch, (exp, cGLC))
#         end
#     end

#     @threads for thid in 1:nthreads()
#         for (exp, cGLC) in Ch

            # ## -------------------------------------------------------------------
            # # handle cache
            # datfile = dat_file(string(DAT_FILE_PREFFIX, method); exp)
            # check_cache(datfile, method) || continue
            
#             ## -------------------------------------------------------------------
#             # SetUp
#             model =  load_model("max_model")
#             biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
#             exp_growth = Fd.val("D", exp)
#             biom_lb, biom_ub = ChU.bounds(model, iJR.BIOMASS_IDER)
#             if biom_ub < exp_growth
#                 lock(WLOCK) do
#                     INDEX[method, :DFILE, exp] = :unfeasible
#                     @info("Not feasible (skipping)", 
#                         biom_ub ,exp_growth, 
#                         thid
#                     ); println()
#                 end
#                 continue
#             end

#             cgD_X = -Fd.cval(:GLC, exp) * Fd.val(:D, exp) / Fd.val(:X, exp)
#             exglc_L = ChU.lb(model, iJR.GLC_EX_IDER)
#             exglc_qta = abs(exglc_L * 0.005)
#             expβ = 0.0
#             vg_avPME = 0.0
#             epouts = Dict()
#             epout = nothing

#             for movround in 1:1000

#                 empty!(epouts)

#                 ## -------------------------------------------------------------------
#                 # GRAD DESCEND
#                 x0 = expβ
#                 x1 = 10.0
#                 maxΔx = max(expβ * 0.05, 1e3)
#                 gdth = 1e-3
#                 target = exp_growth
#                 beta_vec = zeros(size(model, 2))
        
#                 upfrec_time = 50 # secunds
#                 last_uptime = time()
#                 gdit = 1
        
#                 ## -------------------------------------------------------------------
#                 function upfun(beta)
        
#                     beta_vec[biomidx] = beta
#                     epouts[beta] = epout = ChEP.maxent_ep(model; 
#                         beta_vec,
#                         alpha = Inf,
#                         maxiter = 5000,  
#                         epsconv = 1e-4, 
#                         verbose = false, 
#                         solution = epout
#                     )

#                     biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
#                     vg_avPME = ChU.av(model, epout, iJR.GLC_EX_IDER)
        
#                     update = gdit == 1 || abs(last_uptime - time()) > upfrec_time || 
#                         epout.status != ChEP.CONVERGED_STATUS
        
#                     update && lock(WLOCK) do
#                         diff = abs.(exp_growth - biom_avPME)
#                         @info(
#                             "Grad descent... ", 
#                             exp, gdit, 
#                             epout.status, epout.iter, 
#                             biom_avPME, exp_growth, diff, 
#                             (biom_lb, biom_ub),
#                             beta, thid
#                         ); println()
#                         last_uptime = time()
#                     end
                    
#                     gdit += 1
#                     return biom_avPME
#                 end
        
#                 ## -------------------------------------------------------------------
#                 # FIND BETA
#                 expβ = UJL.grad_desc(upfun; x0, x1, th = gdth, maxΔx, 
#                     target, 
#                     maxiter = 5000, 
#                     verbose = false
#                 )
        
#                 ## -------------------------------------------------------------------
#                 # MOVE V_UB
#                 Δstep = 0.5
#                 exglc_lb, exglc_ub = ChU.bounds(model, iJR.GLC_EX_IDER)

#                 # lb is the uptake limit
#                 dist = cgD_X - vg_avPME
#                 Δexglc_lb = sign(dist) * max(exglc_qta, abs(dist * Δstep))
#                 exglc_lb = min(exglc_ub,
#                     min(cgD_X, 
#                         max(exglc_L, exglc_lb + Δexglc_lb)
#                     )
#                 )
#                 ChU.lb!(model, iJR.GLC_EX_IDER, exglc_lb)
                
#                 ## -------------------------------------------------------------------
#                 # INFO AND CONV
#                 biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
#                 gd_err = abs(exp_growth - biom_avPME) / exp_growth
#                 conv = cgD_X <= vg_avPME && epout.status == ChEP.CONVERGED_STATUS && gd_err < gdth
                
#                 lock(WLOCK) do
#                     @info("Round Done", 
#                         movround, conv,
#                         dist, exglc_qta, Δexglc_lb,
#                         (vg_avPME, cgD_X), 
#                         exglc_ub, exglc_lb,  exglc_L, 
#                         thid
#                     ); println()
#                 end
#                 conv && break
        
#             end #  for movround in 1:1000

#             ## -------------------------------------------------------------------
#             lock(WLOCK) do

#                 # Storing
#                 dat = Dict()
#                 dat[:exp_beta] = expβ
#                 dat[:epouts] = epouts
#                 dat[:model] = model |> ChU.compressed_model

#                 # caching
#                 serialize(datfile, dat)
#                 INDEX[method, :DFILE, exp] = datfile

#                 biom_avPME = ChU.av(epouts[expβ])[biomidx]
#                 diff = abs.(exp_growth - biom_avPME)
#                 @info("Finished ",
#                     exp, expβ, 
#                     length(epouts),
#                     biom_avPME, exp_growth, diff, 
#                     thid
#                 ); println()
#             end

#         end # for (exp, cGLC) in Ch
#     end # for thid in 1:nthreads()
# end

# ## -------------------------------------------------------------------
# # Further convergence
# let
#     method = ME_Z_EXPECTED_G_MOVING

#     iterator = Fd.val("cGLC") |> enumerate |> collect
#     @threads for (exp, cGLC) in iterator

#         datfile = INDEX[method, :DFILE, exp]
#         datfile == :unfeasible && continue
#         dat = deserialize(datfile)
#         model, epouts = ChU.uncompressed_model(dat[:model]) , dat[:epouts]

#         exp_growth = Fd.val(:D, exp)
#         exp_beta = dat[:exp_beta]
#         exp_epout = epouts[exp_beta]

#         lock(WLOCK) do
#             @info("Converging...", 
#                 exp, method,
#                 exp_beta, exp_epout.status, 
#                 threadid()
#             ); println()
#         end
#         converg_status = get!(dat, :converg_status, :undone)
#         converg_status == :done && continue
        
#         model = ChLP.fva_preprocess(model; verbose = false, 
#             check_obj = iJR.BIOMASS_IDER
#         )
        
#         new_epout = nothing
#         try
#             biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
#             beta_vec = zeros(size(model, 2)); 
#             beta_vec[biomidx] = exp_beta
#             new_epout = ChEP.maxent_ep(model; 
#                 beta_vec, alpha = Inf, 
#                 epsconv = 1e-5, verbose = false, 
#                 solution = exp_epout, maxiter = 5000
#             )
#         catch err; @warn("ERROR", err); println() end
        
#         biom_avPME = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, iJR.BIOMASS_IDER)
#         fail = isnan(biom_avPME) || biom_avPME == 0.0 
#         epouts[exp_beta] = fail ? exp_epout : new_epout
        
#         # Saving
#         lock(WLOCK) do
#             @info("Saving...", 
#                 exp, method, 
#                 exp_beta, 
#                 new_epout.status,
#                 new_epout.iter,
#                 threadid()
#             ); println()
#         end
#         dat[:model] = ChU.compressed_model(model)
#         dat[:converg_status] = :done
#         serialize(datfile, dat)
#     end
# end

# ## ----------------------------------------------------------------------------
# # ME_Z_EXPECTED_G_BOUNDED
# # initial approach
# let
#     # global setup
#     method = ME_Z_EXPECTED_G_BOUNDED

#     # orig model
#     iterator = Fd.val(:D) |> enumerate |> collect 
#     @threads for (exp, D) in iterator
#         thid = threadid()

        # ## -------------------------------------------------------------------
        # # handle cache
        # datfile = dat_file(string(DAT_FILE_PREFFIX, method); exp)
        # check_cache(datfile, method) || continue

#         # prepare model
#         model = load_model(exp)
#         biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
#         M, N = size(model)
#         exp_growth = Fd.val(:D, exp)
#         growth_ub = ChU.ub(model, iJR.BIOMASS_IDER)
#         feasible = exp_growth < growth_ub
#         biom_lb, biom_ub = ChU.bounds(model, iJR.BIOMASS_IDER)
#         if biom_ub < exp_growth
#             lock(WLOCK) do
#                 INDEX[method, :DFILE, exp] = :unfeasible
#                 @info("Not feasible (skipping)", 
#                     exp, method, 
#                     biom_ub ,exp_growth, 
#                     thid
#                 ); println()
#             end
#             continue
#         end
#         ChU.ub!(model, iJR.BIOMASS_IDER, growth_ub * 1.1) # open a beat helps EP

#         lock(WLOCK) do
#             nzabs_range = ChU.nzabs_range(model.S)
#             @info("Starting... ", 
#                 exp, method,
#                 size(model), nzabs_range, 
#                 feasible,
#                 threadid()
#             ); println()
#         end
#         !feasible && continue

#         # simulation
#         dat = isfile(datfile) ? deserialize(datfile) : Dict()
#         epouts = get!(dat, :epouts, Dict())
#         init_len = length(epouts)
#         beta_vec = zeros(N)
#         approach_status = get!(dat, :approach_status, :running)
#         if approach_status == :finished 
#             lock(WLOCK) do
#                 INDEX[method, :DFILE, exp] = datfile
#             end
#             continue
#         end
#         convth = 0.05
        
#         # log approach
#         epout_seed = isempty(epouts) ? nothing : epouts[maximum(keys(epouts))]
#         betas = [0.0; 10.0.^(3:0.05:15)]
#         nan_beta = first(betas)

#         for approach in [:log_approach, :linear_approach]
            
#             lock(WLOCK) do
#                 @info("Starting", 
#                     exp, method, approach, 
#                     length(epouts),
#                     threadid()
#                 ); println()
#             end

#             for beta in betas

#                 nan_beta = beta
#                 haskey(epouts, beta) && continue

#                 beta_vec[biomidx] = beta
#                 epout = nothing
#                 try
#                     epout = ChEP.maxent_ep(model; 
#                         beta_vec, alpha = Inf, damp = 0.9, epsconv = 1e-4, 
#                         maxvar = 1e50, minvar = 1e-50, verbose = false, solution = epout_seed,
#                         maxiter = 1000
#                     )
#                 catch err; end

#                 # info
#                 biom_avPME = isnothing(epout) ? 0.0 : ChU.av(model, epout, biomidx)
#                 lock(WLOCK) do
#                     @info("Results", exp, method, beta, 
#                         exp_growth, growth_ub, biom_avPME, 
#                         length(epouts),
#                         threadid()
#                     ); println()
#                 end

#                 # error conditions
#                 fail = isnothing(epout) || isnan(biom_avPME) || biom_avPME == 0.0 
#                 fail && break

#                 # updating
#                 epout_seed = epouts[beta] = epout

#                 # convergence
#                 converr = abs(biom_avPME - exp_growth)/exp_growth
#                 conv = converr < convth || biom_avPME > exp_growth 
#                 conv && break

#             end # for betas

#             # Catching
#             update = init_len != length(epouts)
#             update && lock(WLOCK) do
#                 serialize(datfile, dat)
#                 @info("Catching", exp, method,  
#                     length(epouts),
#                     basename(datfile),
#                     threadid()
#                 ); println()
#             end

#             # lineal approach
#             last_beta = maximum(keys(epouts))
#             step = abs(last_beta - nan_beta) / 100.0
#             iszero(step) && break
#             betas = range(last_beta, 1.0e15; step)

#         end # for approach

#         # saving
#         lock(WLOCK) do
#             INDEX[method, :DFILE, exp] = datfile
#             dat[:approach_status] = :finished
#             dat[:model] = model |> ChU.compressed_model
#             dat[:exp_beta] = maximum(keys(epouts))
#             serialize(datfile, dat)
#             @info("Finished", exp, method,  
#                 length(epouts),
#                 threadid()
#             ); println()
#         end
       
#     end # for (exp, D)
# end


# ## ----------------------------------------------------------------------------
# # ME_Z_EXPECTED_G_BOUNDED
# # Further convergence
# let
#     method = ME_Z_EXPECTED_G_BOUNDED

#     iterator = Fd.val(:D) |> enumerate |> collect 
#     @threads for (exp, D) in iterator

#         datfile = INDEX[method, :DFILE, exp]
#         dat = deserialize(datfile)
#         model, epouts = ChU.uncompressed_model(dat[:model]) , dat[:epouts]
        
#         exp_growth = Fd.val(:D, exp)
#         exp_beta = maximum(keys(epouts))
#         exp_epout = epouts[exp_beta]

#         lock(WLOCK) do
#             @info("Converging...", exp, method,
#                 exp_beta, exp_epout.status, 
#                 threadid()
#             ); println()
#         end
#         converg_status = get!(dat, :converg_status, :undone)
#         converg_status == :done && continue

#         new_epout = nothing
#         if exp_epout.status == :unconverged
#             try;
#                 biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
#                 beta_vec = zeros(size(model, 2)); 
#                 beta_vec[biomidx] = exp_beta
#                 new_epout = ChEP.maxent_ep(model; 
#                     beta_vec, alpha = Inf, damp = 0.98, epsconv = 1e-4, 
#                     maxvar = 1e50, minvar = 1e-50, verbose = false, 
#                     solution = exp_epout, maxiter = 5000
#                 )
#             catch err; @warn("ERROR", err); println() end

#             biom_avPME = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, iJR.BIOMASS_IDER)
#             fail = isnan(biom_avPME) || biom_avPME == 0.0 
#             epouts[exp_beta] = fail ? exp_epout : new_epout
#         end
        
#         # Saving
#         lock(WLOCK) do
#             @info("Saving...", exp, exp_beta, 
#                 exp_epout.status,
#                 threadid()
#             ); println()
#         end
#         dat[:converg_status] = :done
#         serialize(datfile, dat)
#     end
# end

# ## ----------------------------------------------------------------------------
# # ME_Z_FIXXED_G_BOUNDED
# let
#     method = ME_Z_FIXXED_G_BOUNDED
#     objider = iJR.BIOMASS_IDER
#     costider = iJR.COST_IDER
#     biomass_f = 0.01

#     iterator = Fd.val(:D) |> enumerate |> collect 
#     @threads for (exp, D) in iterator

#         # handle cache
#         datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
#         if isfile(datfile)
#             lock(WLOCK) do
#                 INDEX[method, :DFILE, exp] = datfile
#                 @info("Cached loaded (skipping)",
#                     exp, D, datfile, threadid()
#                 ); println()
#             end
#             continue
#         end

#         # setup
#         model = load_model(exp)
#         biomidx = ChU.rxnindex(model, objider)
#         M, N = size(model)
#         exp_growth = Fd.val("D", exp)
#         fbaout = ChLP.fba(model, objider, costider)
#         fba_growth = ChU.av(model, fbaout, objider)
#         ub_growth = min(fba_growth, exp_growth)
#         ChU.ub!(model, objider, ub_growth * (1.0 + biomass_f))
#         ChU.lb!(model, objider, ub_growth * (1.0 - biomass_f))
#         model = ChLP.fva_preprocess(model, 
#             check_obj = objider,
#             verbose = false
#         )

#         lock(WLOCK) do
#             @info("Doing... ", 
#                 exp, method, 
#                 D, threadid()
#             ); println()
#         end

#         # maxent
#         epout = ChEP.maxent_ep(model; 
#             alpha = Inf, damp = 0.985, epsconv = 1e-4, 
#             verbose = false, maxiter = 5000
#         )
            
#         # storing
#         lock(WLOCK) do
#             # Storing
#             dat = Dict()
#             dat[:epouts] = Dict(0.0 => epout)
#             dat[:model] = model |> ChU.compressed_model

#             # caching
#             serialize(datfile, dat)
#             INDEX[method, :DFILE, exp] = datfile

#             @info("Finished ", exp, threadid())
#             println()
#         end
#     end
# end

# ## ----------------------------------------------------------------------------
# # ME_Z_OPEN_G_OPEN
# let
#     method = ME_Z_OPEN_G_OPEN
#     objider = iJR.BIOMASS_IDER

#     iterator = Fd.val(:D) |> enumerate |> collect 
#     @threads for (exp, D) in iterator

#         # handle cache
#         datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
#         if isfile(datfile)
#             lock(WLOCK) do
#                 INDEX[method, :DFILE, exp] = datfile
#                 @info("Cached loaded (skipping)",
#                     exp, D, datfile, threadid()
#                 ); println()
#             end
#             continue
#         end

#         # setup
#         model = load_model(exp)
#         biomidx = ChU.rxnindex(model, objider)

#         lock(WLOCK) do
#             @info("Doing... ", 
#                 exp, method, 
#                 D, threadid()
#             ); println()
#         end

#         # maxent
#         epout = ChEP.maxent_ep(model; 
#             alpha = Inf, damp = 0.9, epsconv = 1e-4, 
#             verbose = false, maxiter = 5000
#         )
            
#         # storing
#         lock(WLOCK) do
#             # Storing
#             dat = Dict()
#             dat[:epouts] = Dict(0.0 => epout)
#             dat[:model] = model |> ChU.compressed_model

#             # caching
#             serialize(datfile, dat)
#             INDEX[method, :DFILE, exp] = datfile

#             @info("Finished ", exp, threadid())
#             println()
#         end
#     end
# end

# ## ----------------------------------------------------------------------------
# # save index
# ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = true)