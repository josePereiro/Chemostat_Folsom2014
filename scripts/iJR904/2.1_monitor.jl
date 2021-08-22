using ProjAssistant
@quickactivate

## ----------------------------------------------------------------------------
@time begin
    using Serialization

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    const iJR = ChF.iJR904

    using Plots
    import GR
    !isinteractive() && GR.inline("png")

    using Serialization
    using Base.Threads

    using SimTools
    const SimT = SimTools

end

## ----------------------------------------------------------------------------
let
    monfile = cachedir(iJR, :ME_MONITOR)
    while !isfile(monfile); @warn("Mon file not found", monfile); sleep(10); end
    mon = SimT.OnDiskMonitor(monfile)
    live_proves = Dict()

    mean_min_val, mean_max_val = -1e2, 1e2
    
    fun(v) = clamp(v, mean_min_val, mean_max_val)

    SimT.watch(mon; wt = 15.0) do ddat

        isempty(ddat) && return
        
        for (exp, tdat) in ddat
            method = get(tdat, :method, "")
            
            # Check activity
            old_lprove = get!(live_proves, exp, -1.0)
            new_lprove = get(tdat, :live_prove, -1.0)
            old_lprove == new_lprove && continue
            @info("Doing", exp, old_lprove, new_lprove)
            live_proves[exp] = new_lprove

            for datk in [:round, :gd]
                kdat = get!(tdat, datk, Dict())
                ps = Plots.Plot[]

                # means
                for (avk, limk) in [
                            (:vg_avPME, :cgD_X), 
                            (:biom_avPME, :exp_growth), 
                        ]

                    avdat = get(kdat, avk, [])
                    p = plot(fun.(avdat);  title = string(avk), 
                        xlabel = "iter", ylabel = string(avk), 
                        lw = 3, label = string(avk)
                    )

                    limdat = get(tdat, limk, 0.0)
                    hline!(p, [fun(limdat)]; label = string(limk), 
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

                sfig(iJR, ps, 
                    @fileid, "maxent_monitor", 
                    (;datk, exp, method), 
                    ".png"
                ); @info("Done", exp, datk)

            end # for datk in

        end # for (exp, tdat)

    end
end