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

using Plots
import GR
GR.inline("png")
using Statistics

## -------------------------------------------------------------------
# Tools
fig_path(name) = joinpath(iJR.MODEL_FIGURES_DIR, name)
dat_path(name) = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, name)
function get_closest_beta(exp, D, obj_ider)

    exp_val = Kd.val(:D, exp)
    model = D["model"]
    sepouts = sort(collect(D["epouts"]); by = first)

    closest_beta = nothing
    closest_diff = Inf
    for (beta, epout) in sepouts
        ep_val = ChU.av(model, epout, obj_ider)
        diff = abs(ep_val - exp_val)
        if closest_diff > diff
            closest_diff = diff
            closest_beta = beta
        end
    end
    closest_beta
end

## -------------------------------------------------------------------
# biomass vs beta
function _plot_biomass_vs_beta(SDATA, fpreffix)

    for (exp, D) in SDATA

        p = plot(;xlabel = "beta", ylabel = "biomass", 
            title = string(iJR.PROJ_IDER, " exp_", exp))
        model = D["model"]
        obj_val = iJR.KAYSER_BIOMASS_IDER

        exp_growth = Kd.val(:D, exp)
        fbaout = D["fbaout"]
        fba_growth = ChU.av(model, fbaout, obj_val)

        for (beta, epout) in D["epouts"]
            ep_growth = ChU.av(model, epout, obj_val)
            scatter!(p, [beta], [ep_growth]; 
                color = :red, label = "", alpha = 0.5)
        end

        hline!(p, [fba_growth]; color = :blue, 
            ls = :dash, lw = 2, label = "")
        hline!(p, [exp_growth]; color = :black, 
            ls = :dot, lw = 2, label = "")

        # saving
        fname = fig_path(string(fpreffix, "_exp_", exp, "_biomass_vs_beta.png"))
        savefig(p, fname)
        @show fname
    end
end

## -------------------------------------------------------------------
# stoi err vs beta
function _plot_stoi_err_vs_beta(SDATA, fpreffix)

    for (exp, D) in SDATA

        model = D["model"]
        M, N = size(model)
        sepouts = sort(collect(D["epouts"]); by = first)

        fbaout = D["fbaout"]
        ex_glc = ChU.av(model, fbaout, iJR.GLC_EX_IDER)

        p = plot(title = string(iJR.PROJ_IDER, " exp_", exp), 
            xlabel = "beta", ylabel = "log10( |stoi_err| / |ex_glc| )")

        f(x) = log10(x)
        mat = zeros(length(sepouts), M)
        EPS = 1e-8
        for (i, (beta, epout)) in enumerate(sepouts)
            mat[i, :] = (abs.(ChU.stoi_err(model, epout)) ./ abs(ex_glc)) .+ EPS
        end
        plot!(p, first.(sepouts), f.(mat); 
            label = "", color = :black, alpha = 0.1
        )
        plot!(p, first.(sepouts), f.(mean.(eachrow(mat))); 
            label = "mean", color = :red, alpha = 0.8, 
            lw = 5, ls = :dash
        )
        # saving
        fname = fig_path(string(fpreffix, "_exp_", exp, "_stoi_err_vs_beta.png"))
        savefig(p, fname)
        @show fname
    end
end

## -------------------------------------------------------------------
# biomass correlation
function _plot_full_biomass_correlation(SDATA, fpreffix)
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    fba_plot = plot(;xlabel = "exp biomass", ylabel = "model biomass", 
        title = string(iJR.PROJ_IDER, " FBA"))
    ep_plot = plot(;xlabel = "exp biomass", ylabel = "model biomass", 
        title = string(iJR.PROJ_IDER, " EP"))

    min_, max_ = Inf, -Inf
    for (exp, D) in SDATA

        model = D["model"]
        M, N = size(model)
        sepouts = sort(collect(D["epouts"]); by = first)
        exp_beta = get_closest_beta(exp, D, obj_ider)
        epout = D["epouts"][exp_beta]
        fbaout = D["fbaout"]
        
        # norm
        fba_val = ChU.av(model, fbaout, obj_ider)
        ep_val = ChU.av(model, epout, obj_ider)
        ep_err = sqrt(ChU.va(model, epout, obj_ider))
        exp_val = Kd.val(:D, exp , nothing)
        isnothing(exp_val) && continue
        @show exp_val fba_val ep_val
        scatter!(fba_plot, [exp_val], [fba_val]; 
            label = "", color = :black, alpha = 0.5, ms = 5)

        ep_err = sqrt(ChU.va(model, epout, obj_ider))
        scatter!(ep_plot, [exp_val], [ep_val]; yerr = [ep_err],
            label = "", color = :black, alpha = 0.5, ms = 5)
        
        min_ = minimum([min_, exp_val, fba_val, ep_val])
        max_ = maximum([max_, exp_val, fba_val, ep_val])
    end

    for plt in [fba_plot, ep_plot]
        plot!(plt, [min_,max_], [min_,max_]; label = "", lw = 2, alpha = 0.4, ls = :dash, color = :black)
    end
    p = plot([fba_plot, ep_plot]...; layout = 2)
    fname = fig_path(string(fpreffix, "_biomass_correlation.png"))
    savefig(p, fname)
    @show fname
    nothing
end

## -------------------------------------------------------------------
# full exchs correlation
function _plot_full_exchs_correlation(SDATA, fpreffix)
    # TODO: package this
    Kd_msd_mets = ["GLC", "CO2", "O2", "AC", "NH4"]
    met_map = iJR.load_mets_map()
    exch_map = iJR.load_exch_met_map()
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    fba_exch_plot = plot(;xlabel = "exp exch / |ex_GLC|", ylabel = "model exch / |ex_GLC|", 
        title = string(iJR.PROJ_IDER, " FBA"))
    ep_exch_plot = plot(;xlabel = "exp exch / |ex_GLC|", ylabel = "model exch / |ex_GLC|", 
        title = string(iJR.PROJ_IDER, " EP"))

    min_, max_ = Inf, -Inf
    for (exp, D) in SDATA

        model = D["model"]
        M, N = size(model)
        sepouts = sort(collect(D["epouts"]); by = first)
        exp_beta = get_closest_beta(exp, D, obj_ider)
        epout = D["epouts"][exp_beta]
        fbaout = D["fbaout"]
        
        # norm
        ex_glc = Kd.uval("GLC", exp, nothing)
        isnothing(ex_glc) && continue
        
        for Kd_met in Kd_msd_mets
            iJR_met = met_map[Kd_met]
            iJR_ex = exch_map[iJR_met]
            fba_val = ChU.av(model, fbaout, iJR_ex)  / abs(ex_glc)
            ep_val = ChU.av(model, epout, iJR_ex)  / abs(ex_glc)
            ep_err = sqrt(ChU.va(model, epout, iJR_ex)) / abs(ex_glc)
            exp_val = -Kd.uval(Kd_met, exp , nothing) / abs(ex_glc)
            isnothing(exp_val) && continue
            
            scatter!(fba_exch_plot, [exp_val], [fba_val]; 
                label = "", color = :black, alpha = 0.5)

            scatter!(ep_exch_plot, [exp_val], [ep_val]; 
                yerr = [ep_err],
                label = "", color = :black, alpha = 0.5)
            
            min_ = minimum([min_, exp_val, fba_val, ep_val])
            max_ = maximum([max_, exp_val, fba_val, ep_val])
        end
    end

    for plt in [fba_exch_plot, ep_exch_plot]
        plot!(plt, [min_,max_], [min_,max_]; label = "", lw = 2, alpha = 0.4, ls = :dash, color = :black)
    end
    p = plot([fba_exch_plot, ep_exch_plot]...; layout = 2)
    fname = fig_path(string(fpreffix, "_exchanges_correlation.png"))
    savefig(p, fname)
    @show fname
    nothing
end

## -------------------------------------------------------------------
# full conc correlation
function _plot_full_conc_correlation(SDATA, fpreffix)
    # TODO: package this
    Kd_msd_mets = ["GLC", "AC", "NH4"]
    met_map = iJR.load_mets_map()
    exch_map = iJR.load_exch_met_map()
    obj_ider = iJR.KAYSER_BIOMASS_IDER

    fba_conc_plot = plot(;xlabel = "exp conc / cGLC", ylabel = "model conc / cGLC", 
        title = string(iJR.PROJ_IDER, " FBA"))
    ep_conc_plot = plot(;xlabel = "exp conc / cGLC", ylabel = "model conc / cGLC", 
        title = string(iJR.PROJ_IDER, " EP"))

    min_, max_ = Inf, -Inf
    for (exp, D) in SDATA

        model = D["model"]
        M, N = size(model)
        sepouts = sort(collect(D["epouts"]); by = first)
        exp_beta = get_closest_beta(exp, D, obj_ider)
        epout = D["epouts"][exp_beta]
        fbaout = D["fbaout"]
        exp_xi = Kd.val("xi", exp, nothing)
        isnothing(exp_xi) && continue
        # norm
        cglc = Kd.cval("GLC", exp, nothing)
        isnothing(cglc) && continue

        for Kd_met in Kd_msd_mets
            iJR_met = met_map[Kd_met]
            iJR_ex = exch_map[iJR_met]
            fba_flx = ChU.av(model, fbaout, iJR_ex)
            ep_flx = ChU.av(model, epout, iJR_ex)
            ep_flxerr = sqrt(ChU.va(model, epout, iJR_ex)) 

            # s = c + u*xi
            exp_c = Kd.cval(Kd_met, exp , 0.0)
            exp_s = Kd.sval(Kd_met, exp , nothing)
            isnothing(exp_s) && continue
            exp_s /= cglc
            fba_s = (exp_c + fba_flx * exp_xi) / cglc
            ep_s = (exp_c + ep_flx * exp_xi)/ cglc
            ep_serr = (ep_flxerr * exp_xi) / cglc 

            scatter!(fba_conc_plot, [exp_s], [fba_s]; 
                label = "", color = :black, alpha = 0.5)

            scatter!(ep_conc_plot, [exp_s], [ep_s]; 
                yerr = [ep_serr],
                label = "", color = :black, alpha = 0.5)
            
            min_ = minimum([min_, exp_s, fba_s, ep_s])
            max_ = maximum([max_, exp_s, fba_s, ep_s])
        end
    end

    for plt in [fba_conc_plot, ep_conc_plot]
        plot!(plt, [min_,max_], [min_,max_]; label = "", lw = 2, alpha = 0.4, ls = :dash, color = :black)
    end
    p = plot([fba_conc_plot, ep_conc_plot]...; layout = 2)
    fname = fig_path(string(fpreffix, "_conc_correlation.png"))
    savefig(p, fname)
    @show fname
    nothing
end


for (fpreffix, filename) in [ ("2", "ep_dat_exp"),
                            ("3", "ep_scaled__dat_exp")]

    ## -------------------------------------------------------------------
    # load data
    DATA = Dict()
    for (exp, _) in Kd.val(:D) |> enumerate
        datfile = dat_path(string(filename, exp, ".bson"))
        !isfile(datfile) && continue
        exp, model, epouts, fbaout = ChU.load_data(datfile)
        @show size(model)
        DATA[exp] = Dict("model" => model, "epouts" => epouts, "fbaout" => fbaout)
    end
    SDATA = sort(collect(DATA); by = first)

    ## -------------------------------------------------------------------
    _plot_biomass_vs_beta(SDATA, fpreffix)
    _plot_stoi_err_vs_beta(SDATA, fpreffix)
    _plot_full_biomass_correlation(SDATA, fpreffix)
    _plot_full_exchs_correlation(SDATA, fpreffix)
    _plot_full_conc_correlation(SDATA, fpreffix)

end