let
    method = ME_Z_EXPECTED_G_EXPECTED
    
    ## -------------------------------------------------------------------
    # Monitor
    mon = UJL.OnDiskMonitor(iJR.MODEL_CACHE_DIR, "monitor.jld2")
    UJL.reset!(mon)

    # Feed jobs
    Ch = Channel(nthreads()) do ch
        cGLCs = Fd.val("cGLC")
        for (exp, cGLC)  in enumerate(cGLCs)
            put!(ch, (exp, cGLC))
        end
    end

    @threads for _ in 1:nthreads()
        thid = threadid()
        for (exp, cGLC) in Ch
            exp == 1 && continue # Test


            ## -------------------------------------------------------------------
            # handle cache
            datfile = dat_file(DAT_FILE_PREFFIX; method, exp)
            check_cache(datfile, method) && continue

            ## -------------------------------------------------------------------
            # SetUp
            model =  load_model("max_model")
            M, N = size(model)
            biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            glcidx = ChU.rxnindex(model, iJR.GLC_EX_IDER)
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
            biom_beta = 0.0
            biom_betas = [biom_beta]
            vg_beta = 0.0
            vg_betas = [vg_beta]
            vg_avPME = 0.0
            biom_avPME = 0.0
            biom_diff = 0.0
            vg_diff = 0.0
            beta_vec = zeros(N)
            epouts = Dict()
            epout = nothing
            epout_vgb0 = nothing

            ## -------------------------------------------------------------------
            # TODO: check E[z] = D and E[vg] <= cgG_X at vg_beta = 0.0

            UJL.record!(mon) do dat
                tdat = get!(dat, exp, Dict())
                tdat[:method] = method
                tdat[:cgD_X] = cgD_X
                tdat[:exp_growth] = exp_growth
            end

            ## -------------------------------------------------------------------
            # GRAD DESCEND
            let

                target = [exp_growth, cgD_X * 0.98]
                x0 = [biom_beta, vg_beta]
                maxΔx = [5e2, 1e3]
                x1 = x0 .+ maxΔx .* 0.01
                gdth = 1e-2
                
                upfrec_time = 15
                last_uptime = time()
                gdit = 1

                ## -------------------------------------------------------------------
                function upfun!(gdmodel)

                    biom_beta, vg_beta = UJL.gd_value(gdmodel)
        
                    # Update maxΔx     
                    gdmodel.maxΔx = [
                        max(5e2, abs(biom_beta) * 0.1),
                        max(1e3, abs(vg_beta) * 0.1),
                    ]
                    
                    beta_vec[biomidx] = biom_beta
                    beta_vec[glcidx] = vg_beta
                
                    # maxent
                    epout = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        maxiter = 5000,  
                        epsconv = 1e-3, 
                        verbose = false, 
                        solution = epout
                    )

                    biom_avPME = ChU.av(model, epout, iJR.BIOMASS_IDER)
                    vg_avPME = ChU.av(model, epout, iJR.GLC_EX_IDER)
                    biom_diff = abs(biom_avPME - exp_growth)
                    vg_diff = abs(vg_avPME - cgD_X)

                    update = gdit == 1 || abs(last_uptime - time()) > upfrec_time || 
                        epout.status != ChEP.CONVERGED_STATUS
                    
                    gderr = gdmodel.ϵi
                    update && lock(WLOCK) do
                        @info(
                            "z grad descent... ", 
                            exp, gdit, gderr,
                            epout.status, epout.iter, 
                            (biom_avPME, exp_growth), biom_diff, 
                            (vg_avPME, cgD_X), vg_diff, 
                            (biom_beta, vg_beta), 
                            thid
                        ); println()
                        last_uptime = time()
                    end

                    UJL.record!(mon) do dat
                        tdat = get!(dat, exp, Dict())
                        gddat = get!(tdat, :gd, Dict())
                        UJL.get!push!(gddat; 
                            vg_beta, biom_beta, 
                            biom_avPME, vg_avPME
                        )
                    end
                    
                    gdit += 1
                    return [biom_avPME, vg_avPME]
                end

                ## -------------------------------------------------------------------
                gdmodel = UJL.grad_desc_vec(upfun!; 
                    x0, x1, gdth, maxΔx, 
                    target, maxiter = 5000, 
                    verbose = false
                )
                biom_beta, vg_beta = UJL.gd_value(gdmodel)
            end

            ## -------------------------------------------------------------------
            # COLLECTING
            push!(biom_betas, biom_beta)
            push!(vg_betas, vg_beta)
            epouts[(biom_beta, vg_beta)] = epout

            ## -------------------------------------------------------------------
            # MONITOR
            UJL.record!(mon) do dat
                tdat = get!(dat, exp, Dict())
                rdat = get!(tdat, :round, Dict())
                UJL.get!push!(rdat; 
                    vg_beta, biom_beta, 
                    biom_avPME, vg_avPME
                )
            end

            ## -------------------------------------------------------------------
            lock(WLOCK) do

                # Storing
                dat = Dict()
                dat[:exp_beta] = (biom_beta, vg_beta)
                dat[:epouts] = epouts
                dat[:model] = model |> ChU.compressed_model

                # caching
                serialize(datfile, dat)
                INDEX[method, :DFILE, exp] = datfile

                @info("Finished ",
                    exp, 
                    epout.status, epout.iter, 
                    (biom_avPME, exp_growth), biom_diff, 
                    (vg_avPME, cgD_X), vg_diff, 
                    (biom_beta, vg_beta), 
                    thid
                ); println()
            end

        end # for (exp, cGLC)
    end # for thid
    UJL.reset!(mon)
end