import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid

    # -------------------------------------------------------------------
    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    const iJR = ChF.iJR904
    const Fd = ChF.FolsomData # experimental data
    const Bd = ChF.BegData    # cost data


    # -------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase

    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import ChemostatPlots
    const ChP = ChemostatPlots
    
    import UtilsJL
    const UJL = UtilsJL

    using Serialization

    # -------------------------------------------------------------------
    using Plots, FileIO
    import GR
    GR.inline("png")

end

## -------------------------------------------------------------------
INDEX = ChU.load_data(iJR.MAXENT_VARIANTS_INDEX_FILE; verbose = false);

## -------------------------------------------------------------------
const HOMO = :HOMO
const BOUNDED = :BOUNDED
const EXPECTED = :EXPECTED
const FBA = :FBA

## -------------------------------------------------------------------
fileid = "2.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))
FLX_IDERS = ["GLC", "PYR", "SUCC", "LAC", "FORM", "AC", "O2", "CO2"]

EXPS = 1:4 # experiments that have both concentration and flx data

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = let
    iders = [FLX_IDERS; "D"]
    colors = Plots.distinguishable_colors(length(iders))
    Dict(ider => color for (ider, color) in zip(iders, colors))
end


method_colors = Dict(
    HOMO => :red,
    BOUNDED => :orange,
    EXPECTED => :blue,
)

## -------------------------------------------------------------------
# Collect
DAT = ChU.DictTree()
let 
    objider = iJR.BIOMASS_IDER
    exch_met_map = iJR.load_exch_met_map()
    Kd_mets_map = iJR.load_mets_map()

    for method in [HOMO, EXPECTED, BOUNDED]
        for exp in EXPS
            
            # !haskey(INDEX, method, :DFILE, exp) && continue
            datfile = INDEX[method, :DFILE, exp]
            dat = deserialize(datfile)
            
            model = dat[:model]
            objidx = ChU.rxnindex(model, objider)
            epouts = dat[:epouts]
            exp_beta = maximum(keys(epouts))
            epout = epouts[exp_beta]
            exp_xi = Fd.val(:xi, exp)

            println()
            @info("Doing", exp, method, length(dat[:epouts]), epout.iter);

            # Biomass
            ep_biom = ChU.av(model, epout, objidx)
            ep_std = sqrt(ChU.va(model, epout, objidx))
            Kd_biom = Fd.val("D", exp)
            Kd_biom_err = Fd.err("D", exp)
            
            # store
            DAT[method, :ep   , :flx, objider, exp] = ep_biom
            DAT[method, :eperr, :flx, objider, exp] = ep_std
            DAT[method, :Kd   , :flx, objider, exp] = Kd_biom
            DAT[method, :Kderr, :flx, objider, exp] = Kd_biom_err
            DAT[:Kd   , :flx, objider, exp] = Kd_biom
            DAT[:Kderr, :flx, objider, exp] = Kd_biom_err
            DAT[method, :fva  , :flx, objider, exp] = ChU.bounds(model, objider)
            
            # fuxes
            for Kd_met in FLX_IDERS

                    model_met = Kd_mets_map[Kd_met]
                    model_exch = exch_met_map[model_met]
                    model_exchi = ChU.rxnindex(model, model_exch)

                    ep_av = ChU.av(model, epout, model_exchi)
                    ep_std = sqrt(ChU.va(model, epout, model_exchi))
                    Kd_flx = Fd.val("u$Kd_met", exp)
                    Kd_err = Fd.err("u$Kd_met", exp)
                    
                    DAT[method, :Kd, :flx, Kd_met, exp] = Kd_flx
                    DAT[method, :Kderr, :flx, Kd_met, exp] = Kd_err
                    DAT[:Kd, :flx, Kd_met, exp] = Kd_flx
                    DAT[:Kderr, :flx, Kd_met, exp] = Kd_err
                    DAT[method, :ep, :flx, Kd_met, exp] = ep_av
                    DAT[method, :eperr, :flx, Kd_met, exp] = ep_std
                    
                    DAT[method, :fva , :flx, Kd_met, exp] = ChU.bounds(model, model_exch)

            end

        end # for exp in EXPS
    
    end # for method

end

## -------------------------------------------------------------------
# beta vs stuff
let
    method = EXPECTED
    cGLC_plt = plot(;xlabel = "cGLC", ylabel = "beta")
    D_plt = plot(;xlabel = "D", ylabel = "beta")
    for exp in EXPS 
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        beta = maximum(keys(dat[:epouts]))

        params = (;label = "", color = exp_colors[exp], 
            alpha = 0.7, ms = 7
        )
        cGLC = Fd.val("cGLC", exp)
        D = Fd.val("D", exp)
        scatter!(cGLC_plt, [cGLC], [beta]; params...)
        scatter!(D_plt, [D], [beta]; params...)
    end
    mysavefig([cGLC_plt, D_plt], "beta_vs_stuff"; method)
end

## -------------------------------------------------------------------
# EP biomass corr
let
    objider = iJR.BIOMASS_IDER
    ps = Plots.Plot[]
    for method in [HOMO, EXPECTED, BOUNDED]
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "model biom", ylabel = "exp biom")
        ep_vals = DAT[method, :ep, :flx, objider, EXPS]
        eperr_vals = DAT[method, :eperr, :flx, objider, EXPS]
        Kd_vals = DAT[method, :Kd, :flx, objider, EXPS]
        Kderr_vals = DAT[method, :Kderr, :flx, objider, EXPS]
        color = [exp_colors[exp] for exp in EXPS]
        m, M = myminmax([Kd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, ep_vals, Kd_vals; 
            xerr = eperr_vals,
            yerr = Kderr_vals,
            label = "", color,
            alpha = 0.7, ms = 7,
            xlim = [m - margin, M + margin],
            ylim = [m - margin, M + margin],
        )
        push!(ps, p)
    end
    layout = (1, length(ps))
    mysavefig(ps, "obj_val_ep_corr"; layout)
end

## -------------------------------------------------------------------
# EXPECTED flux vs beta
let
    objider = iJR.BIOMASS_IDER
    method = EXPECTED
    p = plot(title = iJR.PROJ_IDER, xlabel = "beta", ylabel = "biom")
    for exp in EXPS 
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        model = dat[:model]
        objidx = ChU.rxnindex(model, objider)
        epouts = dat[:epouts]
        exp_beta = maximum(keys(epouts))
        exp_xi = Fd.val("xi", exp)
        scatter!(p, [exp_beta], [Fd.val("D", exp)], ms = 12, color = :white, label = "")

        betas = collect(keys(epouts)) |> sort
        bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
        scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

    end
    mysavefig(p, "obj_val_vs_beta"; method)
end

## -------------------------------------------------------------------
# total correlations
let
    dat_prefix = :flx
    iders = FLX_IDERS
    zoom_lim = [-0.5, 1.5]
    
    tot_ps = Plots.Plot[]
    zoom_ps = Plots.Plot[]
    for method in [HOMO, EXPECTED, BOUNDED]                                       
        ep_vals = DAT[method, :ep, dat_prefix, iders, EXPS]
        ep_errs = DAT[method, :eperr, dat_prefix, iders, EXPS]
        Kd_vals = DAT[method, :Kd, dat_prefix, iders, EXPS]
        Kd_errs = DAT[method, :Kderr, dat_prefix, iders, EXPS]
        
        
        diffsign = sign.(Kd_vals) .* sign.(ep_vals)
        Kd_vals = abs.(Kd_vals) .* diffsign
        ep_vals = abs.(ep_vals) .* diffsign

        color = [ider_colors[ider] for ider in iders, exp in EXPS]
        m, M = myminmax([ep_vals; Kd_vals])


        scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
        # ep corr
        p1 = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
            ylabel = "model signdiff $(dat_prefix)",
            xlabel = "exp signdiff $(dat_prefix)", 
        )
        scatter!(p1, Kd_vals, ep_vals; yerr = ep_errs, xerr = Kd_errs, scatter_params...)
        plot!(p1, [m,M], [m,M]; ls = :dash, color = :black, label = "")
        push!(tot_ps, deepcopy(p1))
        
        p2 = plot!(p1; xlim = zoom_lim, ylim = zoom_lim)
        push!(zoom_ps, deepcopy(p1))

    end

    layout = (1, length(tot_ps))
    pname = string(dat_prefix, "_tot_corr")
    mysavefig(tot_ps, pname; layout)
    
    layout = (1, length(zoom_ps))
    pname = string(dat_prefix, "_zoomed_corr")
    mysavefig(zoom_ps, pname; layout)

end

## -------------------------------------------------------------------
# fva bounds
let
   
    ps = Plots.Plot[]
    for ider = FLX_IDERS
        p = plot(title = ider, xlabel = "replica", ylabel = "flx")
        xticks =  (EXPS, string.(EXPS))
        
        Kd_vals = DAT[:Kd, :flx, ider, EXPS]
        plot!(p, EXPS, Kd_vals; 
            label = "exp", color = :black, alpha = 0.8, lw = 3, xticks)

        for method in [HOMO, EXPECTED, BOUNDED]             
            color = method_colors[method]    
            
            ep_vals = DAT[method, :ep, :flx, ider, EXPS]
            plot!(p, EXPS, ep_vals; 
                label = string(method), color, alpha = 0.5, lw = 5, ls = :dash, xticks)
            
            fva_ranges = DAT[method, :fva, :flx, ider, EXPS]
            plot!(p, EXPS, last.(fva_ranges);  
                label = "", color, alpha = 0.8, ls = :dot, lw = 3, xticks)
            plot!(p, EXPS, first.(fva_ranges); 
                label = "", color, alpha = 0.8, ls = :dot, lw = 3, xticks)
        end
        push!(ps, p)
    end
    pname = string("bound_study")
    mysavefig(ps, pname)
    
end

## -------------------------------------------------------------------
# marginal distributions
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    Kd_mets_map = iJR.load_mets_map()
    exch_met_map = iJR.load_exch_met_map()

    # Iders
    model_iders, Kd_iders = [objider], ["D"]
    for Kd_met in FLX_IDERS
        model_met = Kd_mets_map[Kd_met]
        model_exch = exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Kd_iders, string("u", Kd_met))
    end
    
    for (model_ider, Kd_ider) in zip(model_iders, Kd_iders)
        ps = Plots.Plot[]
        ps_bs = Plots.Plot[]
        for exp in EXPS
            p = plot(title = string(Kd_ider, " exp: ", exp))
            p_bs = plot(title = string(Kd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Kd_av = Fd.val(Kd_ider, exp)
            
            # EP
            for method in [BOUNDED, EXPECTED, HOMO]
                color = method_colors[method]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                # fbaout = dat[:fbaout]
                        
                # ChP.plot_marginal!(p, model, [epout, fbaout], model_exch; legend = false)
                ChP.plot_marginal!(p, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.6, lw = 5)
                
                m = minimum([m, ep_av, Kd_av])
                M = maximum([M, ep_av, Kd_av])
                margin = maximum([margin, 3 * ep_va])

                if method == EXPECTED
                    for (beta, epout) in sort(epouts; by = first)
                        ep_av = ChU.av(model, epout, model_ider)
                        ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                        alpha = 0.15
                        color = method_colors[method]
                        ChP.plot_marginal!(p_bs, model, epout, model_ider; 
                            legend = false, color, alpha, lw = 1)

                        if beta == exp_beta
                            ChP.plot_marginal!(p_bs, model, epout, model_ider; 
                                legend = false, color, 
                                alpha = 1.0, lw = 3
                            )
                            break
                        end
                    end
                    push!(ps_bs, p_bs)
                end

            end
            # Experimental
            vline!(p, [Kd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            vline!(p_bs, [Kd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            
            plot!(p; xlim = [m - margin, M + margin], size)
            plot!(p_bs; xlim = [m - margin, M + margin], size)
            push!(ps, p)
        end

        for k in [:xi, :D, :uGLC]
            p = plot(;title = Kd_ider, size)
            xticks =  (EXPS, string.(EXPS))
            vals = [Fd.val(k, exp) for exp in EXPS]
            p = bar!(p, EXPS, vals; title = k, label = "", xticks)
            push!(ps, p)
            push!(ps_bs, p)
        end

        pname = string(Kd_ider, "_marginals")
        mysavefig(ps, pname)

        method = EXPECTED
        pname = string(Kd_ider, "_marginals_vs_beta")
        mysavefig(ps_bs, pname; method)
    end

end 

## -------------------------------------------------------------------
# marginals v2
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    Kd_mets_map = iJR.load_mets_map()
    exch_met_map = iJR.load_exch_met_map()

    # Iders
    model_iders, Kd_iders = [objider], ["D"]
    for Kd_met in FLX_IDERS
        model_met = Kd_mets_map[Kd_met]
        model_exch = exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Kd_iders, string("u", Kd_met))
    end
    
    for (model_ider, Kd_ider) in zip(model_iders, Kd_iders)
        marg_params = (;xlabel = string(Kd_ider), yaxis = nothing, ylabel = "prob")

        epps = Plots.Plot[]
        exps = Plots.Plot[]
        for method in [BOUNDED, EXPECTED, HOMO]
            expp = plot(;title = string("Experimental"), marg_params...)
            epp = plot(;title = string(" MaxEnt: ", method), marg_params...)
            margin, m, M = -Inf, Inf, -Inf
            
            # EP
            for exp in EXPS
                Kd_av = Fd.val(Kd_ider, exp)
                color = exp_colors[exp]    

                datfile = INDEX[method, :DFILE, exp]
                dat = deserialize(datfile)
                model = dat[:model]
                objidx = ChU.rxnindex(model, objider)
                epouts = dat[:epouts]
                exp_beta = maximum(keys(epouts))
                epout = epouts[exp_beta]
                ep_av = ChU.av(model, epout, model_ider)
                ep_va = sqrt(ChU.va(model, epout, model_ider))
                        
                ChP.plot_marginal!(epp, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.8, lw = 3)
                
                m = minimum([m, ep_av, Kd_av])
                M = maximum([M, ep_av, Kd_av])
                margin = maximum([margin, 3 * ep_va])

                # Experimental
                vline!(expp, [Kd_av]; label = "", lw = 3, color, alpha = 0.8)
                
            end
            
            map([expp, epp]) do p
                plot!(p; xlim = [m - margin, M + margin], size)
            end

            push!(epps, epp)
            push!(exps, expp)
        end

        extras = Plots.Plot[]
        for k in [:xi, :D, Kd_ider] |> unique
            p = plot(;title = "Experimental", size, 
                xlabel = "rep", ylabel = string(k))
            xticks =  (EXPS, string.(EXPS))
            vals = [Fd.val(k, exp) for exp in EXPS]
            color = [exp_colors[exp] for exp in EXPS]
            p = bar!(p, EXPS, vals; label = "", xticks, color)
            push!(extras, p)
        end

        ps = Plots.Plot[exps; epps; extras]
        layout = (3, 3)
        pname = string(Kd_ider, "_marginals_v2")
        mysavefig(ps, pname; layout)

    end # for (model_ider, Kd_ider)

end 

## -------------------------------------------------------------------
# leyends
# TODO fix this...
let
    for (title, colors) in [
            ("exp", exp_colors), 
            ("iders", ider_colors),
            ("method", method_colors)
        ]
    p = plot(; framestyle = :none)
        scolors = sort(collect(colors); by = (p) -> string(first(p)))
        for (id, color) in scolors
            scatter!(p, [0], [0];
                thickness_scaling = 1,
                color, ms = 8, label = string(id),
                legendfontsize=10, 
                # size = [300, 900],
                # legend = :left
            )
        end
        mysavefig(p, "$(title)_color_legend")
    end
end