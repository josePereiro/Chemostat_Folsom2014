import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

# Run add https://github.com/josePereiro/Chemostat_Heerden2013.jl in the Julia Pkg REPL to install the
# package, then you must activate the package enviroment (see README)
import Chemostat_Heerden2013
const Hd  = Chemostat_Heerden2013.HeerdenData;
const ECC = Chemostat_Heerden2013.EColiCore

# run add "https://github.com/josePereiro/Chemostat" in the 
# julia Pkg REPL for installing the package
import Chemostat
import Chemostat.LP.MathProgBase
const ChU = Chemostat.Utils
const ChSS = Chemostat.SteadyState
const ChLP = Chemostat.LP
const ChEP = Chemostat.MaxEntEP
const ChSU = Chemostat.SimulationUtils

using Plots

## -------------------------------------------------------------------
# save results
DATA = ChU.load_data(ECC.MAXENT_RES_FILE)
SDATA = sort(collect(DATA); by = first);

## -------------------------------------------------------------------
# flux vs beta
let
    p = plot(title = "EColiCore", xlabel = "beta", ylabel = "biom")
    for (exp, D) in SDATA
        model = D["model"]
        objidx = ChU.rxnindex(model, ECC.BIOMASS_IDER)
        epouts = D["epouts"]
        exp_beta = D["exp_beta"]
        exp_xi = D["exp_xi"]
        scatter!(p, [exp_beta], [Hd.val("D", exp)], ms = 12, color = :white, label = "")

        betas = collect(keys(epouts)) |> sort
        bioms = [ChU.av(model, epouts[beta], objidx) for beta in betas]
        scatter!(p, betas, bioms, label = "", color = :black, alpha = 0.2)

    end
    p
end

## -------------------------------------------------------------------
# beta_exp corr
let
    p = plot(title = "EColiCore", xlabel = "exp biom", ylabel = "model biom")
    for (exp, D) in SDATA

        model = D["model"]
        objidx = ChU.rxnindex(model, ECC.BIOMASS_IDER)
        exp_beta = D["exp_beta"]
        epout = D["epouts"][exp_beta]

        scatter!(p, [ChU.av(model, epout, objidx)], [Hd.val("D", exp)], 
            color = :white, label = "")

    end
    p
end

## -------------------------------------------------------------------
# total correlations
let
    ep_p = plot(title = "EColiCore (EP)", 
        xlabel = "model flux", ylabel = "exp flux")
    fba_p = plot(title = "EColiCore (FBA)", 
        xlabel = "model flux", ylabel = "exp flux")

    m, M = Inf, -Inf
    for (exp, D) in SDATA

        model = D["model"]
        objidx = ChU.rxnindex(model, ECC.BIOMASS_IDER)
        exp_beta = D["exp_beta"]
        epout = D["epouts"][exp_beta]
        fbaout = D["fbaout"]

        ep_xs, ep_errs = [], []
        fba_xs = []
        ys = []

        # Biomass
        push!(fba_xs, ChU.av(model, fbaout, objidx))
        push!(ep_xs, ChU.av(model, epout, objidx))
        push!(ys, Hd.val("D", exp))
        
        # mets
        for Hd_met in Hd.msd_mets
            try
                model_met = ECC.Hd_mets_map[Hd_met]
                model_exch = ECC.exch_met_map[model_met]

                push!(fba_xs, ChU.av(model, fbaout, model_exch))
                push!(ep_xs, ChU.av(model, epout, model_exch))
                push!(ep_errs, sqrt(ChU.va(model, epout, model_exch)))
                push!(ys, Hd.val("u$Hd_met", exp))
            catch 
                @warn string(Hd_met, " fails")
            end
        end

        scatter!(ep_p, ep_xs, ys; xerr = ep_errs, color = :black, label = "")
        scatter!(fba_p, fba_xs, ys, color = :black, label = "")
        m = min(m, minimum(ep_xs), minimum(fba_xs), minimum(ys))
        M = max(M, maximum(ep_xs), maximum(fba_xs), maximum(ys))

    end

    plot!(ep_p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    plot!(fba_p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    plot([ep_p, fba_p]..., size = [700, 400])
end

## -------------------------------------------------------------------
# total correlations (conc)
let
    ep_p = plot(title = "EColiCore (EP)", 
        xlabel = "model conc", ylabel = "exp conc")
    fba_p = plot(title = "EColiCore (FBA)", 
        xlabel = "model conc", ylabel = "exp conc")

    m, M = Inf, -Inf
    for (exp, D) in SDATA

        model = D["model"]
        objidx = ChU.rxnindex(model, ECC.BIOMASS_IDER)
        exp_beta = D["exp_beta"]
        exp_xi = Hd.val("xi", exp)
        epout = D["epouts"][exp_beta]
        fbaout = D["fbaout"]

        ep_xs, ep_errs = [], []
        fba_xs = []
        ys = []
        
        # mets
        for Hd_met in Hd.msd_mets
            try
                model_met = ECC.Hd_mets_map[Hd_met]
                model_exch = ECC.exch_met_map[model_met]

                # conc (s = c + u*xi)
                c = Hd.val("c$Hd_met", exp, 0.0)

                u = ChU.av(model, fbaout, model_exch)
                push!(fba_xs, max(c + u*exp_xi, 0.0))

                u = ChU.av(model, epout, model_exch)
                push!(ep_xs, max(c + u*exp_xi, 0.0))

                uerr = sqrt(ChU.va(model, epout, model_exch)) * exp_xi
                push!(ep_errs, uerr)
                
                push!(ys, Hd.val("s$Hd_met", exp))
                
            catch 
                @warn string(Hd_met, " fails")
            end
        end

        @assert length(ep_xs) == length(fba_xs) == length(ys)
        scatter!(ep_p, ep_xs, ys; xerr = ep_errs, color = :black, label = "")
        scatter!(fba_p, fba_xs, ys, color = :black, label = "")
        m = min(m, minimum(ep_xs), minimum(fba_xs), minimum(ys))
        M = max(M, maximum(ep_xs), maximum(fba_xs), maximum(ys))

    end

    # plot!(ep_p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    # plot!(fba_p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    plot([ep_p, fba_p]..., size = [700, 400])
end