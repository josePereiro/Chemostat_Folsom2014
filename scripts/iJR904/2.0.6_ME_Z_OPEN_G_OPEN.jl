let
    method = ME_Z_OPEN_G_OPEN
    objider = iJR.BIOMASS_IDER

    iterator = Fd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        # handle cache
        datfile = dat_file(DAT_FILE_PREFFIX; exp, method)
        if isfile(datfile)
            lock(WLOCK) do
                INDEX[method, :DFILE, exp] = datfile
                @info("Cached loaded (skipping)",
                    exp, D, datfile, threadid()
                ); println()
            end
            continue
        end

        # setup
        model = load_model(exp)
        biomidx = ChU.rxnindex(model, objider)

        lock(WLOCK) do
            @info("Doing... ", 
                exp, method, 
                D, threadid()
            ); println()
        end

        # maxent
        epout = ChEP.maxent_ep(model; 
            alpha = Inf, damp = 0.9, epsconv = 1e-4, 
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
            INDEX[method, :DFILE, exp] = datfile

            @info("Finished ", exp, threadid())
            println()
        end
    end
end