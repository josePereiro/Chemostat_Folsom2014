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
    method = ME_MAX_POL

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
                model_exch = Fd_rxns_map[ider]
                model_exchi = ChU.rxnindex(model, model_exch)

                model_flx = ChU.av(model, epout, model_exchi)
                exp_flx = Fd.uval(ider, exp)

                sum += (model_flx/glc_flx - exp_flx/glc_flx)^2
                N += 1
            end
            
            push!(MSEs, sum / N)
        end

        scatter!(p, first.(betas), MSEs; color = :black,
            label = "", m = 8, alpha = 0.8
        )
        plot!(p, first.(betas), MSEs; color = :black,
            label = "", ls = :dash, alpha = 0.8
        )
        vline!(p, [first(exp_beta)]; color = :black, 
            label = "", ls = :dot, lw = 3, alpha = 0.9
        )
        push!(ps, p)
    end
    mysavefig(ps, "MSE_vs_beta")
end
