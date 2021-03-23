import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Heerden2013")

@time begin
    ## -------------------------------------------------------------------
    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    const Fd  = ChF.FolsomData;
    const ECC = ChF.EColiCore

    ## -------------------------------------------------------------------
    import Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP
    const ChEP = Chemostat.MaxEntEP
    const ChSU = Chemostat.SimulationUtils
    
    ## -------------------------------------------------------------------
    import UtilsJL
    const UJL = UtilsJL

    using Base.Threads
end

## ----------------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_MAX_POL                = :ME_MAX_POL                 # 
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
# ME_MAX_POL
let
    method = ME_MAX_POL
    
    ## -------------------------------------------------------------------
    # Monitor
    mon = UJL.OnDiskMonitor(ECC.MODEL_CACHE_DIR, "monitor.jld2")
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
            exp != 4 && continue # Test
        end
    end

end

## -------------------------------------------------------------------
# save results
# ChU.save_data(ECC.MAXENT_RES_FILE, DATA)