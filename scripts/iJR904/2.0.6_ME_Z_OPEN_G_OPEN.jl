let
    method = ME_Z_OPEN_G_OPEN
    objider = iJR.BIOMASS_IDER

    iterator = Fd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        # handle cache
        datfile = dat_file(;exp, method)
        # check_cache(;exp, method) && continue

        # setup
        model = iJR.load_model("fva_models", exp)

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
            dat[:exp_beta] = 0.0
            dat[:epouts] = Dict(0.0 => epout)
            dat[:model] = model |> ChU.compressed_model

            # caching
            serialize(datfile, dat)

            @info("Finished ", exp, threadid())
            println()
        end
    end
end