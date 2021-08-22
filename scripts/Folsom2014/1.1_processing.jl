using ProjAssistant
@quickactivate 

# ------------------------------------------------------------------
@time begin

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    const Fd = ChF.FolsomData # experimental data

    using Plots
    import GR
    !isinteractive() && GR.inline("png")
end

## -------------------------------------------------------------------
# D vs X
let
    p = plot(; title = "Kayser", 
        xlabel = string("X (", Fd.unit(:X), ")"), 
        ylabel = string("D (", Fd.unit(:D), ")")
    )
    plot!(p, Fd.val(:D), Fd.val(:X); 
        label = "", ls = :dash, lw = 3, color = :black, 
        alpha = 0.6
    )
    scatter!(p, Fd.val(:D), Fd.val(:X); 
        label = "", m = 8, color = :black
    )
    sfig(Fd, p, 
        @fileid, "X_vs_D", ".png"
    ) 
end

## -------------------------------------------------------------------
# BALANCE
let
    ps = Plots.Plot[]
    for met in [:GLC]
        p = plot(; title = string("Balance: ", met), 
            xlabel = "feed", 
            ylabel = "exch + drain" 
        )

        exps = 1:4

        feed = Fd.cval.(met, exps) .* Fd.val.(:D, exps)
        exch = Fd.uval.(met, exps) .* Fd.val.(:X, exps) .|> abs
        drain = zeros(length(exch)) # This experiments hase sg == 0

        scatter!(p, feed, exch .+ drain; 
            label = "", m = 8, color = :black
        )
        vals = [feed; exch .+ drain] |> sort
        plot!(p, vals, vals;
            label = "", ls = :dash, lw = 3, alpha = 0.6, color = :black
        )
        push!(ps, p)
    end
    sfig(Fd, ps, 
        @fileid, "Balances", ".png"
    ) 
end