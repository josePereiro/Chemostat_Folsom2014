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
            exp_growth = Fd.val("D", exp) # experimental biom growth rate (equals D)
            # glc per biomass unit supply
            cgD_X = -Fd.cval(:GLC, exp) * Fd.val(:D, exp) / Fd.val(:X, exp)
            biom_beta = 0.0 # current biomass beta
            biom_betas = [biom_beta] # biomass beta round history
            vg_beta = 0.0 # current vg beta
            vg_betas = [vg_beta] # vg beta round history
            vg_avPME = 0.0 # current vg exchange average
            vg_avPME_vgb0 = 0.0 # vg exchange average at vg_beta = 0
            biom_avPME = 0.0 # current biom exchange average
            biom_avPME_vgb0 = 0.0 # biom exchange average at vg_beta = 0
            biom_diff = 0.0 # distance between exp_growth and biom_avPME
            vg_diff = 0.0 # distance between cgD_X and vg_diff
            beta_vec = zeros(N) # ep beta vector
            epouts = Dict() # epout pool for each round (will contain the solution)
            epout = nothing # current epout (used as seed)
            hasvalid_moments = false # a flag that indicate is the momentous are valid
            isbeta_stationary = false # a flag that indicates if betas reach stability
            roundconv = false # global (round) converge flag

            epmaxiter = 2000 # maxiter for each maxent_ep
            gdmaxiter = 3000 # maxiter for each gradient descent
            gdth = 0.01  # th of each gradient descend
            roundth = 0.05 # th of the whole simulation
            stw = 10 # beta stability check window
            stth = 0.1 # beta stability check th
            smooth = 0.1 # gd smooth th

            # After a while without converge, accelerate
            turbo_iter0 = 100 # iter for initing turbo
            turbo_frec = 10 # iter frec for apply turbo
            turbo_factor = 2.0 # turbo strength 

            # damp will be reduced when detected
            damp_factor = 0.5 # gd damp penalty factor
            biom_gddamp = 1.0 # biom gd current damp
            vg_gddamp = 1.0 # vg gd current damp

            gdit = -1 # current gd iter
            gderr = -1 # current gd error
            last_uptime = -1 # time to check if gd needs to update
            upfrec_time = 15 # update info frequency

            # the whole simulation converge as ~log, 
            # so I force betas increment by stepping
            beta_scale_len0 = 3 # starting round for beta stepping
            beta_scale_factor = 1.0 # stepping scaling factor

            rounditer = 1 # current round iter
            maxrounds = 50 # max no of rounds

            # monitor
            UJL.record!(mon) do dat
                tdat = get!(dat, exp, Dict())
                tdat[:cgD_X] = cgD_X
                tdat[:method] = method
                tdat[:exp_growth] = exp_growth
            end

            ## -------------------------------------------------------------------
            # closure functions
            function check_roundconv()
                hasvalid_biom_moment = abs(biom_avPME - exp_growth)/abs(exp_growth) <= roundth
                hasvalid_vg_moment = abs(vg_avPME) <= abs(cgD_X) || 
                    abs(vg_avPME - cgD_X)/abs(cgD_X) <= roundth
                hasvalid_moments = hasvalid_biom_moment && hasvalid_vg_moment
                isbeta_stationary = UJL.is_stationary(biom_betas, stth, stw) && 
                    UJL.is_stationary(vg_betas, stth, stw)
                return hasvalid_moments || isbeta_stationary
            end

            ## -------------------------------------------------------------------
            function print_info(msg, mac = @info; varargs...)
                mac(msg, 
                    varargs...,
                    epout.status, epout.iter, 
                    (biom_avPME_vgb0, biom_avPME, exp_growth), biom_diff, 
                    (vg_avPME_vgb0, vg_avPME, cgD_X), vg_diff, 
                    (biom_beta, vg_beta), 
                    thid
                ); println()
            end

            ## -------------------------------------------------------------------
            function gd_core_fun(gdmodel; msg)

                beta_vec[biomidx] = biom_beta
                beta_vec[glcidx] = vg_beta
                
                # MAXENT
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
                    print_info(msg;
                        exp, rounditer, gdit, gderr, 
                        vg_gddamp, biom_gddamp,
                    )
                    last_uptime = time()
                end

                # MONITOR
                UJL.record!(mon) do dat
                    tdat = get!(dat, exp, Dict())
                    tdat[:live_prove] = rand()
                    gddat = get!(tdat, :gd, Dict())
                    UJL.get!push!(gddat; 
                        vg_beta, biom_beta, 
                        biom_avPME, vg_avPME
                    )
                end

                # RUN OUT OF PATIENT
                (gdit >= turbo_iter0 && rem(gdit, turbo_frec) == 0) && (gdmodel.maxΔx *= turbo_factor)
                
                gdit += 1
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
                        gd_core_fun(gdmodel; msg = "z grad descent... ")
                        return biom_avPME
                    end

                    ## -------------------------------------------------------------------
                    function z_break_cond(gdmodel)
                        roundconv = check_roundconv()
                        zconv = abs(biom_avPME - exp_growth)/abs(exp_growth) <= gdth
                        roundconv || zconv
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
                # AT VG BETA 0 MOMENTS
                firstround = rounditer == 1
                firstround && let
                    biom_avPME_vgb0 = biom_avPME
                    vg_avPME_vgb0 = vg_avPME
                end

                ## -------------------------------------------------------------------
                # FORCE VG BOUNDARY
                let
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
                        gd_core_fun(gdmodel; msg = "vg grad descent... ")
                        return vg_avPME
                    end

                    ## -------------------------------------------------------------------
                    function vg_break_cond(epmodel)
                        vgconv = abs(vg_avPME) <= abs(cgD_X)
                        roundconv = check_roundconv()
                        roundconv || vgconv
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
                end

                ## -------------------------------------------------------------------
                # COLLECTING
                push!(biom_betas, biom_beta)
                push!(vg_betas, vg_beta)
                empty!(epouts) # Test
                epouts[(biom_beta, vg_beta)] = epout

                ## -------------------------------------------------------------------
                # CEHCK NAN
                if any(isnan.(biom_betas)) || any(isnan.(vg_betas))
                    print_info("Nan detected (BIG PROBLEMS HERE)", @warn; 
                        exp, rounditer
                    )
                    break
                end
                
                ## -------------------------------------------------------------------
                # MONITOR
                UJL.record!(mon) do dat
                    tdat = get!(dat, exp, Dict())
                    tdat[:live_prove] = rand()
                    rdat = get!(tdat, :round, Dict())
                    UJL.get!push!(rdat; 
                        vg_beta, biom_beta, 
                        biom_avPME, vg_avPME
                    )
                end

                ## -------------------------------------------------------------------
                # PRINT INFO
                lock(WLOCK) do
                    print_info("Round Done"; exp, rounditer, 
                        hasvalid_moments, isbeta_stationary, roundconv
                    )
                end
                
                ## -------------------------------------------------------------------
                # BETA SCALING
                scalebeta = length(biom_betas) >= beta_scale_len0 && 
                    length(vg_betas) >= beta_scale_len0 
                scalebeta && let
                    biom_beta_step = biom_betas[end] - biom_betas[end - 1]
                    biom_beta += biom_beta_step * beta_scale_factor

                    vg_beta_step = vg_betas[end] - vg_betas[end - 1]
                    vg_beta += vg_beta_step * beta_scale_factor
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

                print_info("Finished "; exp, rounditer)
            end

        end # for exp, cGLC
    end # for thid
    UJL.reset!(mon)

end
