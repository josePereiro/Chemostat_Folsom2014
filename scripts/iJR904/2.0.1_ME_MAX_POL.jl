let
    method = ME_MAX_POL
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

            ## -------------------------------------------------------------------
            # handle cache
            datfile = dat_file(DAT_FILE_PREFFIX; method, exp)
            check_cache(datfile, exp, method) && continue

            ## -------------------------------------------------------------------
            # SetUp
            model =  load_model("max_model")
            M, N = size(model)
            biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            glcidx = ChU.rxnindex(model, iJR.GLC_EX_IDER)
            exp_growth = Fd.val("D", exp)

            cgD_X = -Fd.cval(:GLC, exp) * Fd.val(:D, exp) / Fd.val(:X, exp)
            biom_beta = 0.0
            biom_betas = [biom_beta]
            vg_beta = 0.0
            vg_betas = [vg_beta]
            vg_avPME = 0.0
            vg_avPME_vgb0 = 0.0
            biom_avPME = 0.0
            biom_avPME_vgb0 = 0.0
            biom_diff = 0.0
            vg_diff = 0.0
            beta_vec = zeros(N)
            epouts = Dict()
            epout = nothing
            epout_vgb0 = nothing
            hasvalid_epout_moments = false
            isbeta_stationary = false
            roundconv = false

            epmaxiter = 2000
            gdmaxiter = 3000
            gdth = 0.01  # th of each gradient descend
            roundth = 0.05 # th of the whole simulation
            upfrec_time = 15
            stw = 10
            stth = 0.1
            smooth = 0.1

            # This will be reduced every times damping is detected
            damp_factor = 0.5
            biom_gddamp = 1.0
            vg_gddamp = 1.0

            beta_step_len0 = 3
            beta_step_scalef = 1.0

            UJL.record!(mon) do dat
                tdat = get!(dat, exp, Dict())
                tdat[:cgD_X] = cgD_X
                tdat[:method] = method
                tdat[:exp_growth] = exp_growth
            end

            ## -------------------------------------------------------------------
            rounditer = 1
            maxrounds = 50
            function check_roundconv()
                hasvalid_epout_moments = abs(vg_avPME - cgD_X)/abs(cgD_X) <= roundth && 
                    abs(biom_avPME - exp_growth)/abs(exp_growth) <= roundth
                isbeta_stationary = UJL.is_stationary(biom_betas, stth, stw) && 
                    UJL.is_stationary(vg_betas, stth, stw)
                return hasvalid_epout_moments || isbeta_stationary
            end

            while true
                ## -------------------------------------------------------------------
                # Z GRAD DESCEND: Match biomass momentums
                let
                    target = exp_growth
                    x0 = biom_beta
                    maxΔx = max(abs(biom_beta) * 0.05, 5e3)
                    minΔx = maxΔx * 0.01
                    x1 = x0 + maxΔx * 0.01
                    senses = [] # To detect damping
                    check_damp_frec = 10
                    dampth = 0.8
                    maxΔx_reduce_factor = 0.9
                    
                    last_uptime = time()
                    gdit = 1

                    ## -------------------------------------------------------------------
                    function z_fun(gdmodel)

                        biom_beta = UJL.gd_value(gdmodel)
            
                        beta_vec[biomidx] = biom_beta
                        beta_vec[glcidx] = vg_beta
                    
                        # maxent
                        epout = ChEP.maxent_ep(model; 
                            beta_vec,
                            alpha = Inf,
                            maxiter = epmaxiter,  
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
                                exp, rounditer, gdit, gderr,
                                epout.status, epout.iter, 
                                (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                                (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                                (biom_beta, vg_beta), 
                                biom_gddamp,
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
                        return biom_avPME
                    end

                    ## -------------------------------------------------------------------
                    function z_break_cond(gdmodel)
                        roundconv = check_roundconv()
                        zconv = abs(biom_avPME - exp_growth)/abs(exp_growth) <= gdth
                        (rounditer > 10 && roundconv) || zconv
                    end

                    ## -------------------------------------------------------------------
                    gdmodel = UJL.grad_desc(z_fun; 
                        x0, x1, gdth, minΔx, maxΔx, smooth,
                        target, maxiter = gdmaxiter, 
                        damp_factor, damp = biom_gddamp,
                        break_cond = z_break_cond,
                        verbose = false
                    )
                    biom_beta = UJL.gd_value(gdmodel)
                    biom_gddamp = gdmodel.damp
                end

                ## -------------------------------------------------------------------
                # CHECK VG VALIDITY
                firstround = rounditer == 1
                firstround && let
                    beta_vec[biomidx] = biom_beta
                    beta_vec[glcidx] = 0.0

                    # maxent
                    epout_vgb0 = ChEP.maxent_ep(model; 
                        beta_vec,
                        alpha = Inf,
                        maxiter = epmaxiter,  
                        epsconv = 1e-3, 
                        verbose = false,
                        solution = epout
                    )   
                    biom_avPME_vgb0 = ChU.av(model, epout_vgb0, biomidx)
                    vg_avPME_vgb0 = ChU.av(model, epout_vgb0, glcidx)
                end

                ## -------------------------------------------------------------------
                # Force vg boundary
                let
                    if abs(vg_avPME_vgb0) <= abs(cgD_X)
                        vg_beta = 0.0
                        epout = epout_vgb0
                        biom_avPME = biom_avPME_vgb0
                        vg_avPME = vg_avPME_vgb0
                        biom_diff = abs(biom_avPME - exp_growth)
                        vg_diff = abs(vg_avPME - cgD_X)
                    else
                        ## -------------------------------------------------------------------
                        # VG GRAD DESCEND: Match biomass momentums
                        target = cgD_X * 0.99 # force to be inside
                        x0 = vg_beta
                        maxΔx = max(abs(vg_beta) * 0.05, 1e3)
                        minΔx = maxΔx * 0.01
                        x1 = x0 + maxΔx * 0.01
                
                        last_uptime = time()
                        gdit = 1
        
                        ## -------------------------------------------------------------------
                        function vg_fun(gdmodel)
        
                            vg_beta = UJL.gd_value(gdmodel)
                
                            beta_vec[biomidx] = biom_beta
                            beta_vec[glcidx] = vg_beta
                        
                            epout = ChEP.maxent_ep(model; 
                                beta_vec,
                                alpha = Inf,
                                maxiter = epmaxiter,  
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
                                    "vg grad descent... ", 
                                    exp, rounditer, gdit, gderr,
                                    epout.status, epout.iter, 
                                    (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                                    (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                                    (biom_beta, vg_beta), 
                                    vg_gddamp,
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
                            return vg_avPME
                        end

                        ## -------------------------------------------------------------------
                        function vg_break_cond(epmodel)
                            vgconv = abs(vg_avPME) <= abs(cgD_X)
                            roundconv = check_roundconv()
                            (rounditer > 10 && roundconv) || vgconv
                        end

                        ## -------------------------------------------------------------------
                        gdmodel = UJL.grad_desc(vg_fun; 
                            x0, x1, gdth, minΔx, maxΔx,
                            break_cond = vg_break_cond,
                            damp_factor, damp = vg_gddamp,
                            target, maxiter = gdmaxiter, 
                            verbose = false
                        )
                        vg_beta = UJL.gd_value(gdmodel)
                        vg_gddamp = gdmodel.damp

                    end # if abs(vg_avPME_vgb0) <= abs(cgD_X)
                end

                ## -------------------------------------------------------------------
                # COLLECTING
                push!(biom_betas, biom_beta)
                push!(vg_betas, vg_beta)
                empty!(epouts) # Test
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
                # PRINTED INFO
                lock(WLOCK) do
                    @info("Round Done", 
                        exp, rounditer, 
                        isbeta_stationary, hasvalid_epout_moments,
                        roundconv,
                        epout.status, epout.iter, 
                        (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                        (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                        (biom_beta, vg_beta), 
                        thid
                    ); println()
                end
                
                ## -------------------------------------------------------------------
                # BETA SCALING
                scalebeta = length(biom_betas) >= beta_step_len0 && 
                    length(vg_betas) >= beta_step_len0 
                scalebeta && let
                    biom_beta_step = biom_betas[end] - biom_betas[end - 1]
                    biom_beta += biom_beta_step * beta_step_scalef

                    vg_beta_step = vg_betas[end] - vg_betas[end - 1]
                    vg_beta += vg_beta_step * beta_step_scalef
                end

                ## -------------------------------------------------------------------

                # BREAK
                roundconv && break
                rounditer += 1
                rounditer > maxrounds && break
            end # round while

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
                    exp, rounditer,  
                    epout.status, epout.iter, 
                    (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                    (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                    (biom_beta, vg_beta), 
                    thid
                ); println()
            end

        end # for exp, cGLC
    end # for thid
    UJL.reset!(mon)

end
