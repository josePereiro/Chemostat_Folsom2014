import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

# ------------------------------------------------------------------
using MAT

import Chemostat
const ChLP = Chemostat.LP
const ChSS = Chemostat.SteadyState
const ChSU = Chemostat.SimulationUtils
const ChU = Chemostat.Utils

import Chemostat_Kayser2005: KayserData, iJR904
const iJR = iJR904
const Kd = KayserData

import ProgressMeter: ProgressThresh, Progress, next!, finish!, update!
using Plots

## ------------------------------------------------------------------
# Biomass medium sensibility
let
    model = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
    obj_ider = iJR.KAYSER_BIOMASS_IDER
    xi = Kd.val("xi") |> minimum
    intake_info = iJR.load_base_intake_info()
    results = Dict()
    factors = 0.0:0.01:1.0
    prog = Progress(length(factors) * length(intake_info))
    for (exch, dat) in intake_info
        c0 = dat["c"]
        res = get!(results, exch, [])
        for f in factors
            dat["c"] = f * c0
            ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)
            growth = try 
                    fbaout = ChLP.fba(model, obj_ider, iJR.COST_IDER)
                    ChU.av(model, fbaout, obj_ider)
                catch err; 0.0 end
            push!(res, growth)

            next!(prog; showvalues = [
                    (:exch, exch),
                    (:f, f),
                    (:c0, c0),
                    (:c, dat["c"]),
                    (:growth, growth)
                ]
            )
        end 
        dat["c"] = c0 
    end
    finish!(prog)

    p = plot(
            title = "Biomass medium sensivility", 
            xlabel = "fraction of initial conc", 
            ylabel = "growth"
        )
    sresults = sort(collect(results); by = (x) -> sum(x[2]))
    lcount = 4
    for (exch, res) in sresults
        lb_ = lcount > 0 ? exch : ""
        plot!(p, factors, res; label = lb_, lw = 3)
        lcount -= 1
    end
    p
end

## ------------------------------------------------------------------
# Checking fba_obj_val < exp_obj_val
let 
    to_map = Kd.val("D") |> enumerate
    for (Di, D) in to_map

        model = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)

        ## Open intakes except Glucose
        intake_info = iJR.load_base_intake_info()
        intake_info[iJR.GLC_EX_IDER]["c"] = Kd.val(:cGLC, Di)
        
        # impose constraint
        xi = Kd.val(:xi, Di)
        ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

        fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER, iJR.COST_IDER);
        fba_obj_val = ChU.av(model, fbaout, iJR.KAYSER_BIOMASS_IDER)
        fba_obj_val = ChU.av(model, fbaout, iJR.KAYSER_BIOMASS_IDER)
        fba_ex_glc_val = ChU.av(model, fbaout, iJR.GLC_EX_IDER)
        fba_ex_glc_b = ChU.bounds(model, iJR.GLC_EX_IDER)
        exp_obj_val = Kd.val("D", Di)

        ChU.tagprintln_inmw("FBA SOLUTION", 
            "\nxi:                      ", xi,
            "\nobj_ider:                ", iJR.KAYSER_BIOMASS_IDER,
            "\nfba fba_ex_glc_val:      ", fba_ex_glc_val,
            "\nfba fba_ex_glc_b:        ", fba_ex_glc_b,
            "\nfba obj_val:             ", fba_obj_val,
            "\nexp obj_val:             ", exp_obj_val,
            "\ncost_ider:               ", iJR.COST_IDER,
            "\nfba cost_val:            ", ChU.av(model, fbaout, iJR.COST_IDER),
            "\n\n"
        )
        (fba_obj_val < exp_obj_val) && @warn "fba objval < exp objval" fba_obj_val exp_obj_val
    end
end

## ------------------------------------------------------------------
# Find cGLC that fit experimental growth
let
    to_map = Kd.val("D") |> enumerate
    for (Di, D) in to_map
        Kd_cGLC = Kd.val(:cGLC, Di)
        Kd_growth = Kd.val(:D, Di)
        xi = Kd.val(:xi, Di)
        model = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
        intake_info = iJR.load_base_intake_info()
        for (exch, info) in intake_info
            info["c"] = iJR.MAX_CONC # open medium
        end
        
        function work_fun(cGLC)
            ## Open intakes except Glucose
            intake_info[iJR.GLC_EX_IDER]["c"] = first(cGLC)
            
            # impose constraint
            ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

            ## fba
            fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER, iJR.COST_IDER)
            fba_growth = ChU.av(model, fbaout, iJR.KAYSER_BIOMASS_IDER)
            return [fba_growth]
        end

        cGLC = ChSU.grad_desc(work_fun; x0 = [Kd_cGLC], x1 = [Kd_cGLC * 0.9], th = 1e-5, 
            C = [Kd_cGLC * 0.1], target = [Kd_growth], maxiters = 500) |> first
        @info "Results" Di Kd_cGLC cGLC
    end

end

## ------------------------------------------------------------------
# Testing scaled model
let
    model = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
    fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER, iJR.COST_IDER);
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", iJR.KAYSER_BIOMASS_IDER,
        "\nsize:             ", size(model),
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJR.KAYSER_BIOMASS_IDER),
        "\ncost_ider:        ", iJR.COST_IDER,
        "\nfba cost_val:     ", ChU.av(model, fbaout, iJR.COST_IDER),
        "\n\n"
    )
    model = ChU.well_scaled_model(model, 100.0; verbose = false)
    fbaout = ChLP.fba(model, iJR.KAYSER_BIOMASS_IDER, iJR.COST_IDER);
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", iJR.KAYSER_BIOMASS_IDER,
        "\nsize:             ", size(model),
        "\nfba obj_val:      ", ChU.av(model, fbaout, iJR.KAYSER_BIOMASS_IDER),
        "\ncost_ider:        ", iJR.COST_IDER,
        "\nfba cost_val:     ", ChU.av(model, fbaout, iJR.COST_IDER),
        "\n\n"
    )
end

## ------------------------------------------------------------------
# old vs new biomass
let

    model0 = ChU.load_data(iJR.BASE_MODEL_FILE; verbose = false)
    to_map = Kd.val("D") |> enumerate
    p = plot(;xlabel = "xi", ylabel = "biomass")

    dat = Dict()
    obj_ider = iJR.KAYSER_BIOMASS_IDER
    for (Di, D) in to_map

        model = deepcopy(model0)
        obj_idx = ChU.rxnindex(model, obj_ider)

        ## Open intakes except Glucose
        intake_info = iJR.load_base_intake_info()
        intake_info[iJR.GLC_EX_IDER]["c"] = Kd.val(:cGLC, Di)
        
        # impose constraint
        xi = Kd.val(:xi, Di)
        ChSS.apply_bound!(model, xi, intake_info; 
            emptyfirst = true, ignore_miss = true)

        fbaout = ChLP.fba(model, obj_ider);
        fba_obj_val = ChU.av(model, fbaout, obj_ider)
        fba_avs = get!(dat, obj_ider, [])
        push!(fba_avs, fba_obj_val)
        exp_obj_val = Kd.val("D", Di)

    end
    plot!(p, Kd.val(:xi), dat[obj_ider]; label = obj_ider, lw = 3, alpha = 0.5)
    scatter!(p, Kd.val(:xi), Kd.val(:D); color = :black, label = "exp")
    p
end
