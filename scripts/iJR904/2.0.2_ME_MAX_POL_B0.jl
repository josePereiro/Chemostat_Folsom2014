function do_max_pol_b0(method, model_key)

    # ep params
    me_params = lglob(iJR, :maxent, :params)
    @extract me_params: alpha epsconv maxiter 
    @extract me_params: damp maxvar minvar

    # handle cache
    is_cached(;method) && return
    
    # setup
    model = iJR.load_model(model_key)

    lock(WLOCK) do
        @info("Doing... ", 
            method, threadid()
        ); println()
    end

    # maxent
    epout = ChEP.maxent_ep(model; 
        alpha, damp, maxiter,
        epsconv, maxvar, minvar, 
    )
        
    # storing
    lock(WLOCK) do
        
        # Storing
        dat = Dict()
        dat[:exp_beta] = 0.0
        dat[:epouts] = Dict(0.0 => epout)
        dat[:model] = model |> ChU.compressed_model

        # caching
        datfile = dat_file(;method)
        sdat(iJR, dat, datfile)

        @info("Finished ", method); println()
    end
end

