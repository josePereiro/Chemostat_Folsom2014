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
        Fderr_vals = DAT[method, :Fderr, :flx, objider, EXPS]
        color = [exp_colors[exp] for exp in EXPS]
        m, M = myminmax([Fd_vals; ep_vals])
        margin = abs(M - m) * 0.1
        scatter!(p, Fd_vals, ep_vals; 
            xerr = Fderr_vals,
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
