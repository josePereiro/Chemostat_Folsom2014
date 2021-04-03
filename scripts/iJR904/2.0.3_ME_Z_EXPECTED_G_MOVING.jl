let
    method = ME_Z_EXPECTED_G_MOVING

    # Feed jobs
    Ch = Channel(1) do ch
        cGLCs = Fd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    @threads for thid in 1:nthreads()
        for (exp, cGLC) in Ch

            ## -------------------------------------------------------------------
            # handle cache
            datfile = dat_file(string(DAT_FILE_PREFFIX, method); exp)
            check_cache(datfile, exp, method) || continue
            
            ## -------------------------------------------------------------------
            # SetUp
            model =  load_model("max_model")
            biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            exp_growth = Fd.val("D", exp)
            biom_lb, biom_ub = ChU.bounds(model, iJR.BIOMASS_IDER)
            if biom_ub < exp_growth
                lock(WLOCK) do
                    INDEX[method, :DFILE, exp] = :unfeasible
                    @info("Not feasible (skipping)", 
                        biom_ub ,exp_growth, 
                        thid
                    ); println()
                end
                continue
            end

            cgD_X = -Fd.cval(:GLC, exp) * Fd.val(:D, exp) / Fd.val(:X, exp)
            exglc_L = ChU.lb(model, iJR.EX_GLC_IDER)
            exglc_qta = abs(exglc_L * 0.005)
            expβ = 0.0
            vg_avPME = 0.0
            epouts = Dict()
            epout = nothing

            for movround in 1:1000

                empty!(epouts)

                ## -------------------------------------------------------------------
                # GRAD DESCEND
                x0 = expβ
                x1 = 10.0
                maxΔx = max(expβ * 0.05, 1e3)
                gdth = 1e-3
                target = exp_growth
                beta_vec = zeros(size(model, 2))
        
                upfrec_time = 50 # secunds
                last_uptime = time()
                gdit = 1
        
                ## -------------------------------------------------------------------
                function upfun(beta)
        
                    beta_vec[biomidx] = beta
                    epouts[beta] = epout = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        maxiter = 5000,  
                        epsconv = 1e-4, 
                        verbose = false, 
                        solution = epout
                    )

                    biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
                    vg_avPME = ChU.av(model, epout, iJR.EX_GLC_IDER)
        
                    update = gdit == 1 || abs(last_uptime - time()) > upfrec_time || 
                        epout.status != ChEP.CONVERGED_STATUS
        
                    update && lock(WLOCK) do
                        diff = abs.(exp_growth - biom_avPME)
                        @info(
                            "Grad descent... ", 
                            exp, gdit, 
                            epout.status, epout.iter, 
                            biom_avPME, exp_growth, diff, 
                            (biom_lb, biom_ub),
                            beta, thid
                        ); println()
                        last_uptime = time()
                    end
                    
                    gdit += 1
                    return biom_avPME
                end
        
                ## -------------------------------------------------------------------
                # FIND BETA
                expβ = UJL.grad_desc(upfun; x0, x1, th = gdth, maxΔx, 
                    target, 
                    maxiter = 5000, 
                    verbose = false
                )
        
                ## -------------------------------------------------------------------
                # MOVE V_UB
                Δstep = 0.5
                exglc_lb, exglc_ub = ChU.bounds(model, iJR.EX_GLC_IDER)

                # lb is the uptake limit
                dist = cgD_X - vg_avPME
                Δexglc_lb = sign(dist) * max(exglc_qta, abs(dist * Δstep))
                exglc_lb = min(exglc_ub,
                    min(cgD_X, 
                        max(exglc_L, exglc_lb + Δexglc_lb)
                    )
                )
                ChU.lb!(model, iJR.EX_GLC_IDER, exglc_lb)
                
                ## -------------------------------------------------------------------
                # INFO AND CONV
                biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
                gderr = abs(exp_growth - biom_avPME) / exp_growth
                conv = cgD_X <= vg_avPME && epout.status == ChEP.CONVERGED_STATUS && gderr < gdth
                
                lock(WLOCK) do
                    @info("Round Done", 
                        movround, conv,
                        dist, exglc_qta, Δexglc_lb,
                        (vg_avPME, cgD_X), 
                        exglc_ub, exglc_lb,  exglc_L, 
                        thid
                    ); println()
                end
                conv && break
        
            end #  for movround in 1:1000

            ## -------------------------------------------------------------------
            lock(WLOCK) do

                # Storing
                dat = Dict()
                dat[:exp_beta] = expβ
                dat[:epouts] = epouts
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)
                INDEX[method, :DFILE, exp] = datfile

                biom_avPME = ChU.av(epouts[expβ])[biomidx]
                diff = abs.(exp_growth - biom_avPME)
                @info("Finished ",
                    exp, expβ, 
                    length(epouts),
                    biom_avPME, exp_growth, diff, 
                    thid
                ); println()
            end

        end # for (exp, cGLC) in Ch
    end # for thid in 1:nthreads()
end

# Further convergence
let
    method = ME_Z_EXPECTED_G_MOVING

    iterator = Fd.val("cGLC") |> enumerate |> collect
    @threads for (exp, cGLC) in iterator

        datfile = INDEX[method, :DFILE, exp]
        datfile == :unfeasible && continue
        dat = deserialize(datfile)
        model, epouts = ChU.uncompressed_model(dat[:model]) , dat[:epouts]

        exp_growth = Fd.val(:D, exp)
        exp_beta = dat[:exp_beta]
        exp_epout = epouts[exp_beta]

        lock(WLOCK) do
            @info("Converging...", 
                exp, method,
                exp_beta, exp_epout.status, 
                threadid()
            ); println()
        end
        converg_status = get!(dat, :converg_status, :undone)
        converg_status == :done && continue
        
        model = ChLP.fva_preprocess(model; verbose = false, 
            check_obj = iJR.BIOMASS_IDER
        )
        
        new_epout = nothing
        try
            biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            beta_vec = zeros(size(model, 2)); 
            beta_vec[biomidx] = exp_beta
            new_epout = ChEP.maxent_ep(model; 
                beta_vec, alpha = Inf, 
                epsconv = 1e-5, verbose = false, 
                solution = exp_epout, maxiter = 5000
            )
        catch err; @warn("ERROR", err); println() end
        
        biom_avPME = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, iJR.BIOMASS_IDER)
        fail = isnan(biom_avPME) || biom_avPME == 0.0 
        epouts[exp_beta] = fail ? exp_epout : new_epout
        
        # Saving
        lock(WLOCK) do
            @info("Saving...", 
                exp, method, 
                exp_beta, 
                new_epout.status,
                new_epout.iter,
                threadid()
            ); println()
        end
        dat[:model] = ChU.compressed_model(model)
        dat[:converg_status] = :done
        serialize(datfile, dat)
    end
end
