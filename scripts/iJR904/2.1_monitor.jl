import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    using Serialization

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    const iJR = ChF.iJR904

    using Plots

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads

end

## ----------------------------------------------------------------------------
fileid = "2.0"
mysavefig(p, pname; params...) = 
    UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
## ----------------------------------------------------------------------------
let
    mon = UJL.OnDiskMonitor(iJR.MODEL_CACHE_DIR, "monitor.jld2")
    UJL.sync_from_disk!(mon)
    UJL.get_cache(mon)
    UJL.watch(mon) do ddat

        # vg_beta, biom_beta, 
        # biom_avPME, vg_avPME
        for (exp, tdat) in ddat
            method = get(tdat, :method, "")
            for datk in [:round, :gd]
                kdat = get!(tdat, datk, Dict())
                ps = Plots.Plot[]

                # means
                for (avk, limk) in [
                            (:vg_avPME, :cgD_X), 
                            (:biom_avPME, :exp_growth), 
                        ]
                    avdat = get(kdat, avk, [])
                    p = plot(avdat;  title = string(avk), 
                        xlabel = "iter", ylabel = string(avk), 
                        lw = 3, label = string(avk)
                    )
                    limdat = get(tdat, limk, 0.0)
                    hline!(p, [limdat]; label = string(limk), 
                        lw = 3, ls = :dash, color = :black, 
                    )
                    plot!(p; legend = :topleft)
                    push!(ps, p)
                end
                
                # betas
                for bk in [:vg_beta, :biom_beta]
                    dat = get(kdat, bk, [])
                    p = plot(dat;  title = string(bk), 
                        xlabel = "iter", ylabel = string(bk), 
                        lw = 3, label = string(bk)
                    )
                    plot!(p; legend = :topleft)
                    push!(ps, p)
                end

                mysavefig(ps, "monitor"; datk, exp, method)
                @info("Done", exp, datk)
            end
        end
    end
end