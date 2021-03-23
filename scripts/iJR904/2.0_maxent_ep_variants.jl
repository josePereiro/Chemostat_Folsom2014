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
    UJL.set_cache_dir(iJR.MODEL_CACHE_DIR)
end

## ----------------------------------------------------------------------------
# globals
const WLOCK = ReentrantLock()
const SIM_GLOBAL_ID = "iJR904_MAXENT_VARIANTS"
const DAT_FILE_PREFFIX =  "maxent_ep_dat_"

const INDEX = UJL.DictTree()
function dat_file(name; kwargs...)
    fname = UJL.mysavename(name, "jls"; kwargs...)
    joinpath(iJR.MODEL_PROCESSED_DATA_DIR, fname)
end

## ----------------------------------------------------------------------------
function load_model(exp::Int, modelkey::String = "fva_models")
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS[modelkey][exp]
    model |> ChU.uncompressed_model
end

function load_model(modelkey::String = "max_model")
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS[modelkey]
    model |> ChU.uncompressed_model
end

## ----------------------------------------------------------------------------
const ME_Z_OPEN_G_OPEN          = :ME_Z_OPEN_G_OPEN           # Do not use extra constraints
const ME_MAX_POL                = :ME_MAX_POL                 # 
const ME_Z_EXPECTED_G_EXPECTED  = :ME_Z_EXPECTED_G_EXPECTED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_BOUNDED   = :ME_Z_EXPECTED_G_BOUNDED    # Match ME and Dy biom average and constraint av_ug
const ME_Z_EXPECTED_G_MOVING    = :ME_Z_EXPECTED_G_MOVING     # 
const ME_Z_FIXXED_G_BOUNDED     = :ME_Z_FIXXED_G_BOUNDED      # Fix biom around observed

## -------------------------------------------------------------------
function check_cache(datfile, exp, method)
    thid = threadid()
    if isfile(datfile)
        lock(WLOCK) do
            INDEX[method, :DFILE, exp] = datfile
            @info("Cached loaded (skipping)",
                exp, datfile, thid
            )
            println()
        end
        return true
    end
    return false
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

# ## ----------------------------------------------------------------------------
# save index
ChU.save_data(iJR.MAXENT_VARIANTS_INDEX_FILE, INDEX; verbose = true)