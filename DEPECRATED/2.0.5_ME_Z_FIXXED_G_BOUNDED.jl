let
    method = ME_Z_FIXXED_G_BOUNDED
    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    biomass_f = 0.01

    iterator = Fd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        check_cache(datfile, exp, method) || continue

        # setup
        model = load_model(exp)
        biomidx = ChU.rxnindex(model, objider)
        M, N = size(model)
        exp_growth = Fd.val("D", exp)
        fbaout = ChLP.fba(model, objider, costider)
        fba_growth = ChU.av(model, fbaout, objider)
        ub_growth = min(fba_growth, exp_growth)
        ChU.ub!(model, objider, ub_growth * (1.0 + biomass_f))
        ChU.lb!(model, objider, ub_growth * (1.0 - biomass_f))
        model = ChLP.fva_preprocess(model, 
            check_obj = objider,
            verbose = false
        )

        lock(WLOCK) do
            @info("Doing... ", 
                exp, method, 
                D, threadid()
            ); println()
        end

        # maxent
        epout = ChEP.maxent_ep(model; 
            alpha = Inf, damp = 0.985, epsconv = 1e-4, 
            verbose = false, maxiter = 5000
        )
            
        # storing
        lock(WLOCK) do
            # Storing
            dat = Dict()
            dat[:epouts] = Dict(0.0 => epout)
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)
            

            @info("Finished ", exp, threadid())
            println()
        end
    end
end