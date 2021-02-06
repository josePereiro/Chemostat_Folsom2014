import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

import SparseArrays
import Base.Threads: @threads, threadid, SpinLock

# -------------------------------------------------------------------
using Serialization

import Chemostat_Kayser2005
import Chemostat_Kayser2005.Chemostat
import Chemostat_Kayser2005.Chemostat.LP: MathProgBase
const ChK = Chemostat_Kayser2005

const iJR = ChK.iJR904
const Kd = ChK.KayserData # experimental data
const Bd = ChK.BegData    # cost data

const Ch = ChK.Chemostat
const ChU = Ch.Utils
const ChSS = Ch.SteadyState
const ChLP = Ch.LP
const ChEP = Ch.MaxEntEP
const ChSU = Ch.SimulationUtils

## -------------------------------------------------------------------

## -------------------------------------------------------------------
let
    # setup
    cache_dir = iJR.MODEL_CACHE_DIR

    # orig model
    for (exp, D) in Kd.val(:D) |> enumerate

        model = ChU.load_data(iJR.BASE_MODEL_FILE)
        model = ChU.fix_dims(model)
        ChU.check_dims(model)
        obj_ider = iJR.KAYSER_BIOMASS_IDER
        fbaout = ChLP.fba(model, obj_ider)
        fba_growth = ChU.av(model, fbaout, obj_ider)
        @info "Orig model" exp size(model) ChU.nzabs_range(model.S) fba_growth
        
        # Scaling
        scale_factor = 1000.0
        model = ChU.well_scaled_model(model, scale_factor)
        fbaout = ChLP.fba(model, obj_ider)
        fba_growth = ChU.av(model, fbaout, obj_ider)
        @info "Scaled model" exp scale_factor size(model) ChU.nzabs_range(model.S) fba_growth

        # prepare model
        exp_growth = Kd.val(:D, exp)
        exp_xi = Kd.val(:xi, exp)
        intake_info = iJR.intake_info(exp)
        ChSS.apply_bound!(model, exp_xi, intake_info; 
            emptyfirst = true)
        fbaout = ChLP.fba(model, obj_ider)
        fba_growth = ChU.av(model, fbaout, obj_ider)
        @info "base model" exp size(model) ChU.nzabs_range(model.S) fba_growth

        # fva model
        fvafile = joinpath(cache_dir, string("fva_cache_exp_", exp, "_sf_", scale_factor, ".bson"))
        if isfile(fvafile)
            model = ChU.load_data(fvafile)
        else
            model = ChLP.fva_preprocess(model, check_obj = obj_ider, verbose = true)
            ChU.save_data(fvafile, model)
        end
        fbaout = ChLP.fba(model, obj_ider)
        fba_growth = ChU.av(model, fbaout, obj_ider)
        @info "fva model" exp size(model) ChU.nzabs_range(model.S) fba_growth

        # simulation
        M, N = size(model)
        objidx = ChU.rxnindex(model, obj_ider)
        datfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, string("ep_scaled__dat_exp", exp, ".bson"))
        epouts = isfile(datfile) ? ChU.load_data(datfile).epouts : Dict()
        beta_vec = zeros(N)

        # log approach
        betas = [0.0; 10.0.^(3:0.05:15)]
        epout_seed = isempty(epouts) ? nothing : epouts[maximum(keys(epouts))]
        nan_beta = first(betas)

        for approach in [:log_approach, :linear_approach]
            
            @info "Starting" exp approach

            for beta in betas

                nan_beta = beta
                haskey(epouts, beta) && continue

                beta_vec[objidx] = beta
                epout = try; 
                    ChEP.maxent_ep(model; 
                            beta_vec, alpha = Inf, damp = 0.9, epsconv = 1e-4, 
                            maxvar = 1e50, minvar = 1e-50, verbose = true, solution = epout_seed,
                            maxiter = 1000
                        )
                catch err; (@warn(err); nothing); end

                isnothing(epout) && break

                # info
                ep_growth = ChU.av(model, epout, objidx)
                @info "Results" exp beta exp_growth fba_growth ep_growth
                println()

                # out conditions
                isnan(ep_growth) && break
                ep_growth == 0.0 && break

                # updating
                epout_seed = epout
                epouts[beta] = epout
            end

            # lineal approach
            last_beta = maximum(keys(epouts))
            bstep = abs(last_beta - nan_beta) / 100.0
            @info "Dev" last_beta nan_beta bstep
            betas = range(last_beta, 1.0e15; step = bstep)
        end

        # fba
        fbaout = let
            lmodel = deepcopy(model)
            ChU.ub!(lmodel, obj_ider, Kd.val(:D, exp))
            ChLP.fba(lmodel, obj_ider, iJR.COST_IDER)
        end

        # saving
        ChU.save_data(datfile, (;exp, model, epouts, fbaout))

    end # for (exp, D)

end