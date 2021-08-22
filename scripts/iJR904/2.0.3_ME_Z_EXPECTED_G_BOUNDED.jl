# initial approach
function do_z_expected_ug_bounded(method, model_key)

    # orig model
    iterator = Fd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator
        thid = threadid()
        
        ## -------------------------------------------------------------------
        # Excluded because of unfeasibility
        exp in IGNORE_EXPS && continue 

        ## -------------------------------------------------------------------
        # handle cache
        is_cached(;method, exp) && continue

        ## -------------------------------------------------------------------
        # ep
        me_params = lglob(iJR, :maxent, :params)
        @extract me_params: alpha epsconv maxiter damp maxvar minvar

        # approach
        convth = 0.05
        converr = nothing

        # prepare model
        model = iJR.load_model(model_key, exp)
        biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        M, N = size(model)
        exp_growth = Fd.val(:D, exp)
        growth_ub = ChU.ub(model, iJR.BIOMASS_IDER)
        ChU.ub!(model, biomidx, growth_ub * 1.1) # open a bit helps EP

        lock(WLOCK) do
            nzabs_range = ChU.nzabs_range(model.S)
            @info("Starting... ", 
                exp, method,
                size(model), 
                nzabs_range, 
                thid
            ); println()
        end

        # simulation globals
        epouts = Dict()
        beta_vec = zeros(N)
        epout_seed = nothing
        epout = nothing
        betas = [0.0; 10.0.^(3:0.05:15)]
        nan_beta = first(betas)

        for approach in [:log_approach, :linear_approach]
            
            lock(WLOCK) do
                @info("Starting", approach, 
                    exp, method,  
                    length(epouts),
                    thid
                ); println()
            end

            for beta in betas

                nan_beta = beta
                haskey(epouts, beta) && continue

                beta_vec[biomidx] = beta
                try
                    epout = ChEP.maxent_ep(model; 
                        beta_vec, alpha, damp, maxiter,
                        epsconv, maxvar, minvar, 
                        verbose = false,
                        solution = epout_seed
                    )
                catch err
                    lock(WLOCK) do
                        @info("Log approach errors", exp, beta, err, thid); 
                        println()
                    end
                    break 
                end
                
                # stop conditions
                biom_avPME = ChU.av(model, epout, biomidx)
                converr = abs(biom_avPME - exp_growth) / exp_growth

                # Check failing
                if (isnan(biom_avPME) || biom_avPME == 0.0)
                    lock(WLOCK) do
                        @info("Log approach fails", exp, beta, biom_avPME, thid); 
                        println()
                    end
                    break
                end

                # info
                lock(WLOCK) do
                    @info("Results", approach, 
                        exp, method, beta, 
                        exp_growth, growth_ub, biom_avPME, 
                        converr, convth,
                        length(epouts),
                        thid
                    ); println()
                end

                # updating
                epout_seed = epouts[beta] = epout

                # convergence
                conv = converr < convth || biom_avPME > exp_growth 
                conv && break

            end # for betas

            # setup lineal approach
            isempty(epouts) && error("LOG FAILS TOO SOON")
            last_beta = maximum(keys(epouts))
            step = abs(last_beta - nan_beta) / 100.0
            iszero(step) && break
            betas = range(last_beta, 1.0e15; step)

        end # for approach

        # further conv
        let
            exp_beta = maximum(keys(epouts))

            @info("Further conv", 
                exp, method,  
                length(epouts),
                exp_beta,
                thid
            ); println()

            beta_vec[biomidx] = exp_beta
            epouts[exp_beta] = ChEP.maxent_ep(model; 
                beta_vec, alpha, damp, 
                maxiter = max(maxiter * 3, 5000),
                epsconv = epsconv * 0.1, 
                maxvar, minvar, 
                verbose = false,
                solution = epouts[exp_beta]
            )
        end

        # saving
        lock(WLOCK) do
            datfile = dat_file(;method, exp)
            dat = Dict()
            dat[:model] = model |> ChU.compressed_model
            dat[:epouts] = epouts
            dat[:exp_beta] = maximum(keys(epouts))
            sdat(iJR, dat, datfile)

            @info("Finished", exp, method,  
                length(epouts),
                thid
            ); println()
        end
       
    end # for (exp, D)

end


