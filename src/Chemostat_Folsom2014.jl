module Chemostat_Folsom2014

    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP

    import UtilsJL
    UtilsJL.gen_top_proj(@__MODULE__)

    include("Utils/Utils.jl")
    include("FolsomData/FolsomData.jl")
    include("BegData/BegData.jl")
    include("EColiCore/EColiCore.jl")
    include("iJR904/iJR904.jl")

    function __init__()
        UtilsJL.create_proj_dirs(@__MODULE__)
    end

end
