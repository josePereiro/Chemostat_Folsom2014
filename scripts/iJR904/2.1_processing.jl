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

# -------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN        = :ME_Z_OPEN_G_OPEN
const ME_Z_EXPECTED_G_MOVING  = :ME_Z_EXPECTED_G_MOVING  
const ME_Z_FIXXED_G_BOUNDED   = :ME_Z_FIXXED_G_BOUNDED
const ME_Z_EXPECTED_G_BOUNDED = :ME_Z_EXPECTED_G_BOUNDED

ALL_METHODS = [ME_Z_OPEN_G_OPEN, ME_Z_FIXXED_G_BOUNDED, 
    ME_Z_EXPECTED_G_BOUNDED, ME_Z_EXPECTED_G_MOVING]

# -------------------------------------------------------------------
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
    ME_Z_OPEN_G_OPEN => :red,
    ME_Z_FIXXED_G_BOUNDED => :orange,
    ME_Z_EXPECTED_G_MOVING => :purple,
    ME_Z_EXPECTED_G_BOUNDED => :blue,
)

exch_met_map = iJR.load_exch_met_map()
Fd_mets_map = iJR.load_mets_map()

## -------------------------------------------------------------------
# Collect
DAT = ChU.DictTree()
let 
    
    DATfile = joinpath(iJR.MODEL_PROCESSED_DATA_DIR, "2.1_DAT.jls")
    # CACHE
    if isfile(DATfile) 
        global DAT = deserialize(DATfile) 
        @info("DAT CACHE LOADED")
        return
    end

    objider = iJR.BIOMASS_IDER
    DAT[:FLX_IDERS] = FLX_IDERS
    DAT[:EXPS] = []

    # Find exps
    for exp in 1:13
        ok = false
        for method in ALL_METHODS
            ok = haskey(INDEX, method, :DFILE, exp) &&
                INDEX[method, :DFILE, exp] != :unfeasible
            !ok && break
        end
        !ok && continue
        push!(DAT[:EXPS], exp)
    end

    for exp in DAT[:EXPS], method in ALL_METHODS
            
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
        Fd_biom = Fd.val("D", exp)
        Fd_biom_err = Fd.err("D", exp)
        
        # store
        DAT[method, :ep   , :flx, objider, exp] = ep_biom
        DAT[method, :eperr, :flx, objider, exp] = ep_std
        DAT[method, :Fd   , :flx, objider, exp] = Fd_biom
        DAT[method, :Fderr, :flx, objider, exp] = Fd_biom_err
        DAT[:Fd   , :flx, objider, exp] = Fd_biom
        DAT[:Fderr, :flx, objider, exp] = Fd_biom_err
        DAT[method, :fva  , :flx, objider, exp] = ChU.bounds(model, objider)
        
        # fluxes
        for Fd_met in FLX_IDERS

                model_met = Fd_mets_map[Fd_met]
                model_exch = exch_met_map[model_met]
                model_exchi = ChU.rxnindex(model, model_exch)

                proj = ChLP.projection2D(model, objider, model_exchi; l = 50)
                ep_av = ChU.av(model, epout, model_exchi)
                ep_std = sqrt(ChU.va(model, epout, model_exchi))
                Fd_flx = Fd.val("u$Fd_met", exp)
                Fd_err = Fd.err("u$Fd_met", exp)
                
                DAT[method, :ep, :proj, Fd_met, exp] = proj
                DAT[method, :Fd, :flx, Fd_met, exp] = Fd_flx
                DAT[method, :Fderr, :flx, Fd_met, exp] = Fd_err
                DAT[:Fd, :flx, Fd_met, exp] = Fd_flx
                DAT[:Fderr, :flx, Fd_met, exp] = Fd_err
                DAT[method, :ep, :flx, Fd_met, exp] = ep_av
                DAT[method, :eperr, :flx, Fd_met, exp] = ep_std
                
                DAT[method, :fva , :flx, Fd_met, exp] = ChU.bounds(model, model_exch)

        end # for Fd_met
    end # for exp in EXPS, for method

    DAT[:EXPS] |> unique! |> sort!
    serialize(DATfile, DAT)
end;
EXPS = DAT[:EXPS]

## -------------------------------------------------------------------
# Inter project comunication
let
    CORR_DAT = isfile(iJR.CORR_DAT_FILE) ? ChU.load_data(iJR.CORR_DAT_FILE) : Dict()
    CORR_DAT[:MAXENT_EP] = DAT
    ChU.save_data(iJR.CORR_DAT_FILE, CORR_DAT)
end

## -------------------------------------------------------------------
# MSE per method
let

    p = plot(;xlabel = "experiment", ylabel = "MSE")
    for method in ALL_METHODS
        MSEs = []
        for exp in EXPS

            sum = 0.0
            N = 0
            glc_flx = DAT[method, :Fd, :flx, "GLC", exp]
            for ider in FLX_IDERS
                model_val = DAT[method, :ep, :flx, ider, exp]
                exp_val = DAT[method, :Fd, :flx, ider, exp]
                sum += (model_val/glc_flx - exp_val/glc_flx)^2
                N += 1
            end
            push!(MSEs, sum / N)
        end
        scatter!(p, EXPS, MSEs; color = method_colors[method],
            label = string(method), m = 8, alpha = 0.8, 
            legend = :topleft
        )
        plot!(p, EXPS, MSEs; color = method_colors[method],
            label = "", ls = :dash, alpha = 0.8
        )
    end
    mysavefig(p, "MSE_per_method")
end

## -------------------------------------------------------------------
# MSE per ider
let
    p = plot(;xlabel = "experiment", ylabel = "MSE")
    for method in ALL_METHODS
        MSEs = []

        for ider in FLX_IDERS
            sum = 0.0
            N = 0
            for exp in EXPS
                glc_flx = DAT[method, :Fd, :flx, "GLC", exp]
                model_val = DAT[method, :ep, :flx, ider, exp]
                exp_val = DAT[method, :Fd, :flx, ider, exp]
                sum += (model_val/glc_flx - exp_val/glc_flx)^2
                N += 1
            end
            push!(MSEs, sum / N)
        end

        scatter!(p, FLX_IDERS, MSEs; color = method_colors[method],
            label = string(method), m = 8, alpha = 0.8, 
            legend = :topleft
        )
        plot!(p, FLX_IDERS, MSEs; color = method_colors[method],
            label = "", ls = :dash, alpha = 0.8
        )
    end
    mysavefig(p, "MSE_per_ider")
end

## -------------------------------------------------------------------
# MSE per beta
let
    method = ME_Z_EXPECTED_G_BOUNDED

    ps = Plots.Plot[]
    for exp in EXPS
        p = plot(;title = string("exp: ", exp), 
            xlabel = "beta", ylabel = "MSE"
        )

        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        epouts = dat[:epouts]
        betas = epouts |> keys |> collect |> sort
        exp_beta = maximum(betas) # dat[:exp_beta]
        model = dat[:model]
        
        MSEs = []
        for beta in betas
            epout = epouts[beta]

            sum = 0.0
            N = 0

            glc_flx = Fd.uval(:GLC, exp)
            for ider in FLX_IDERS

                model_met = Fd_mets_map[ider]
                model_exch = exch_met_map[model_met]
                model_exchi = ChU.rxnindex(model, model_exch)

                model_flx = ChU.av(model, epout, model_exchi)
                exp_flx = Fd.uval(ider, exp)

                sum += (model_flx/glc_flx - exp_flx/glc_flx)^2
                N += 1
            end
            
            push!(MSEs, sum / N)
        end

        scatter!(p, betas, MSEs; color = :black,
            label = "", m = 8, alpha = 0.8
        )
        plot!(p, betas, MSEs; color = :black,
            label = "", ls = :dash, alpha = 0.8
        )
        vline!(p, [exp_beta]; color = :black, 
            label = "", ls = :dot, lw = 3, alpha = 0.9
        )
        push!(ps, p)
    end
    mysavefig(ps, "MSE_vs_beta")
end

## -------------------------------------------------------------------
# proj 2D
let
    method = ME_Z_EXPECTED_G_MOVING
    biom_ider = iJR.BIOMASS_IDER

    ps_pool = Dict()
    for exp in EXPS

        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
            
        model = dat[:model]
        for Fd_ider in FLX_IDERS

            # 2D Projection
            p = plot(;title = string("Folsom2014, exp: ", exp), 
                xlabel = string(biom_ider), ylabel = string(Fd_ider),
                legend = :left
            )
            proj = DAT[method, :ep, :proj, Fd_ider, exp]
            ChP.plot_projection2D!(p, proj; l = 50)

            # cgD/X
            input = -Fd.cval(Fd_ider, exp, 0.0) * Fd.val(:D, exp) / Fd.val(:X, exp)
            hline!(p, [input]; lw = 3, color = :black, ls = :solid, label = "input")

            # EXPERIMENTAL FLXS
            exp_biom = DAT[method, :Fd, :flx, biom_ider, exp]
            exp_biom_err = DAT[method, :Fderr, :flx, biom_ider, exp]
            exp_exch = DAT[method, :Fd, :flx, Fd_ider, exp]
            exp_exch_err = DAT[method, :Fderr, :flx, Fd_ider, exp]
            scatter!(p, [exp_biom], [exp_exch]; 
                xerr = [exp_biom_err], yerr = [exp_exch_err],
                m = 8, color = :red, label = "exp"
            )
            
            # MAXENT FLXS
            ep_biom = DAT[method, :ep, :flx, biom_ider, exp]
            ep_biom_err = DAT[method, :eperr, :flx, biom_ider, exp]
            ep_exch = DAT[method, :ep, :flx, Fd_ider, exp]
            ep_exch_err = DAT[method, :eperr, :flx, Fd_ider, exp]
            scatter!(p, [ep_biom], [ep_exch];
                xerr = [ep_biom_err], yerr = [ep_exch_err],
                m = 8, color = :blue, label = "maxent"
            )

            # mysavefig(p, "polytope"; Fd_ider, exp, method)
            ps_pool[(exp, Fd_ider)] = deepcopy(p)
        end
    end

    # collect 
    for exp in EXPS
        ps = Plots.Plot[ps_pool[(exp, Fd_ider)] for Fd_ider in FLX_IDERS]
        mysavefig(ps, "polytope"; exp, method)
    end

    for Fd_ider in FLX_IDERS
        ps = Plots.Plot[ps_pool[(exp, Fd_ider)] for exp in EXPS]
        mysavefig(ps, "polytope"; Fd_ider, method)
    end
end

## -------------------------------------------------------------------
# beta vs stuff
let
    method = ME_Z_EXPECTED_G_MOVING
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
    for method in ALL_METHODS
        p = plot(title = string(iJR.PROJ_IDER, " method: ", method), 
            xlabel = "exp biom", ylabel = "model biom")
        ep_vals = DAT[method, :ep, :flx, objider, EXPS]
        eperr_vals = DAT[method, :eperr, :flx, objider, EXPS]
        Fd_vals = DAT[method, :Fd, :flx, objider, EXPS]
        Kderr_vals = DAT[method, :Fderr, :flx, objider, EXPS]
        color = [exp_colors[exp] for exp in EXPS]
        m, M = myminmax([Fd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, Fd_vals, ep_vals; 
            xerr = Kderr_vals,
            yerr = eperr_vals,
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
# flux vs beta
let
    objider = iJR.BIOMASS_IDER
    method = ME_Z_EXPECTED_G_MOVING
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
    for method in ALL_METHODS                                       
        ep_vals = DAT[method, :ep, dat_prefix, iders, EXPS]
        ep_errs = DAT[method, :eperr, dat_prefix, iders, EXPS]
        Fd_vals = DAT[method, :Fd, dat_prefix, iders, EXPS]
        Fd_errs = DAT[method, :Fderr, dat_prefix, iders, EXPS]
        
        
        diffsign = sign.(Fd_vals) .* sign.(ep_vals)
        Fd_vals = abs.(Fd_vals) .* diffsign
        ep_vals = abs.(ep_vals) .* diffsign

        color = [ider_colors[ider] for ider in iders, exp in EXPS]
        m, M = myminmax([ep_vals; Fd_vals])


        scatter_params = (;label = "", color, ms = 7, alpha = 0.7)
        # ep corr
        p1 = plot(title = "$(iJR.PROJ_IDER) (EP) $method", 
            xlabel = "exp signdiff $(dat_prefix)", 
            ylabel = "model signdiff $(dat_prefix)",
        )
        scatter!(p1, Fd_vals, ep_vals; yerr = ep_errs, xerr = Fd_errs, scatter_params...)
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
        
        Fd_vals = DAT[:Fd, :flx, ider, EXPS]
        plot!(p, EXPS, Fd_vals; 
            label = "exp", color = :black, alpha = 0.8, lw = 3, xticks)

        for method in [ME_Z_OPEN_G_OPEN, ME_Z_EXPECTED_G_BOUNDED, ME_Z_FIXXED_G_BOUNDED]             
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
    Fd_mets_map = iJR.load_mets_map()
    exch_met_map = iJR.load_exch_met_map()

    # Iders
    model_iders, Fd_iders = [objider], ["D"]
    for Fd_met in FLX_IDERS
        model_met = Fd_mets_map[Fd_met]
        model_exch = exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Fd_iders, string("u", Fd_met))
    end
    
    for (model_ider, Fd_ider) in zip(model_iders, Fd_iders)
        ps = Plots.Plot[]
        ps_bs = Plots.Plot[]
        for exp in EXPS
            p = plot(title = string(Fd_ider, " exp: ", exp))
            p_bs = plot(title = string(Fd_ider, " exp: ", exp))
            margin, m, M = -Inf, Inf, -Inf
            Fd_av = Fd.val(Fd_ider, exp)
            
            # EP
            for method in ALL_METHODS
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
                        
                ChP.plot_marginal!(p, model, [epout], model_ider; 
                    legend = false, color, alpha = 0.6, lw = 5)
                
                m = minimum([m, ep_av, Fd_av])
                M = maximum([M, ep_av, Fd_av])
                margin = maximum([margin, 3 * ep_va])

                if method == ME_Z_EXPECTED_G_MOVING
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
            vline!(p, [Fd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            vline!(p_bs, [Fd_av]; label = "", lw = 6, color = :black, alpha = 0.3)
            
            plot!(p; xlim = [m - margin, M + margin], size)
            plot!(p_bs; xlim = [m - margin, M + margin], size)
            push!(ps, p)
        end

        for k in [:xi, :D, :uGLC]
            p = plot(;title = Fd_ider, size)
            xticks =  (EXPS, string.(EXPS))
            vals = [Fd.val(k, exp) for exp in EXPS]
            p = bar!(p, EXPS, vals; title = k, label = "", xticks)
            push!(ps, p)
            push!(ps_bs, p)
        end

        # legend
        p = plot(;title = "Legend", size)
        for (method, color) in method_colors
            color = method_colors[method]
            bar!(p, [string(method)], [1]; yaxis = nothing, 
                color, xrotation = 35, label = ""
            )
        end
        push!(ps, p)
        push!(ps_bs, p)

        pname = string(Fd_ider, "_marginals")
        mysavefig(ps, pname)

        method = ME_Z_EXPECTED_G_MOVING
        pname = string(Fd_ider, "_marginals_vs_beta")
        mysavefig(ps_bs, pname; method)
    end

end 

## -------------------------------------------------------------------
# marginals v2
let 
    objider = iJR.BIOMASS_IDER
    size = [300, 250]
    Fd_mets_map = iJR.load_mets_map()
    exch_met_map = iJR.load_exch_met_map()

    # Iders
    model_iders, Fd_iders = [objider], ["D"]
    for Fd_met in FLX_IDERS
        model_met = Fd_mets_map[Fd_met]
        model_exch = exch_met_map[model_met]
        push!(model_iders, model_exch)
        push!(Fd_iders, string("u", Fd_met))
    end
    
    for (model_ider, Fd_ider) in zip(model_iders, Fd_iders)
        marg_params = (;xlabel = string(Fd_ider), yaxis = nothing, ylabel = "prob")

        epps = Plots.Plot[]
        exps = Plots.Plot[]
        for method in ALL_METHODS
            expp = plot(;title = string("Experimental"), marg_params...)
            epp = plot(;title = string(" MaxEnt: ", method), marg_params...)
            margin, m, M = -Inf, Inf, -Inf
            
            # EP
            for exp in EXPS
                Fd_av = Fd.val(Fd_ider, exp)
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
                
                m = minimum([m, ep_av, Fd_av])
                M = maximum([M, ep_av, Fd_av])
                margin = maximum([margin, 3 * ep_va])

                # Experimental
                vline!(expp, [Fd_av]; label = "", lw = 3, color, alpha = 0.8)
                
            end
            
            map([expp, epp]) do p
                plot!(p; xlim = [m - margin, M + margin], size)
            end

            push!(epps, epp)
            push!(exps, expp)
        end

        extras = Plots.Plot[]
        for k in [:xi, :X, Fd_ider] |> unique
            p = plot(;title = "Experimental", size, 
                xlabel = "rep", ylabel = string(k))
            xticks =  (EXPS, string.(EXPS))
            vals = [Fd.val(k, exp) for exp in EXPS]
            color = [exp_colors[exp] for exp in EXPS]
            p = bar!(p, EXPS, vals; label = "", xticks, color)
            push!(extras, p)
        end

        # legend
        p = plot(;title = "Legend", size)
        for (method, color) in method_colors
            color = method_colors[method]
            bar!(p, [string(method)], [1]; yaxis = nothing, 
                color, xrotation = 35, label = ""
            )
        end
        push!(extras, p)

        ps = Plots.Plot[exps; epps; extras]
        layout = (3, 4)
        pname = string(Fd_ider, "_marginals_v2")
        mysavefig(ps, pname; layout)

    end # for (model_ider, Fd_ider)

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