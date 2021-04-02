module Chemostat_Folsom2014

    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_top_proj(@__MODULE__)

    include("FolsomData/FolsomData.jl")
    include("Utils/Utils.jl")
    include("BegData/BegData.jl")
    include("iJR904/iJR904.jl")
    # include("EColiCore/EColiCore.jl")

    function __init__()
        UJL.create_proj_dirs(@__MODULE__)
    end

end
