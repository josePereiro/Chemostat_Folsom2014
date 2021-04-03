import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    using Serialization

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    const iJR = ChF.iJR904
    const Fd = ChF.FolsomData # experimental data
    const Bd = ChF.BegData    # cost data

    import Chemostat
    import Chemostat.LP: MathProgBase
    const Ch = ChF.Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    using Statistics

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads
    UJL.set_cache_dir(iJR.cachedir())
end

## ----------------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const DAT_FILE_PREFFIX = "maxent_ep_dat"

## ----------------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_MAX_POL                = :ME_MAX_POL                 # 
const ME_MAX_POL_B0             = :ME_MAX_POL_B0                 # 
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
function dat_file(;kwargs...)
    fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
    iJR.procdir(fname)
end

function check_cache(;kwargs...)
    thid = threadid()
    datfile = dat_file(;kwargs...)
    isfile(datfile) && lock(WLOCK) do
        @info("Cached loaded (skipping)",
            datfile, thid
        ); println()
    end; isfile(datfile)
end

## -------------------------------------------------------------------
# ME_MAX_POL
include("2.0.1_ME_MAX_POL.jl")

# ## -------------------------------------------------------------------
# # ME_Z_EXPECTED_G_EXPECTED
# include("2.0.2_ME_Z_EXPECTED_G_EXPECTED.jl")

# ## -------------------------------------------------------------------
# # # ME_Z_EXPECTED_G_MOVING
# include("2.0.3_ME_Z_EXPECTED_G_MOVING.jl")

# ## ----------------------------------------------------------------------------
# include("2.0.4_ME_Z_EXPECTED_G_BOUNDED.jl")

# ## ----------------------------------------------------------------------------
# # ME_Z_FIXXED_G_BOUNDED
# include("2.0.5_ME_Z_FIXXED_G_BOUNDED.jl")

# ## ----------------------------------------------------------------------------
# # ME_Z_OPEN_G_OPEN
# include("2.0.6_ME_Z_OPEN_G_OPEN.jl")

## ----------------------------------------------------------------------------
# ME_MAX_POL_B0
include("2.0.7_ME_MAX_POL_B0.jl")