module Chemostat_Folsom2014

    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP

    include("Utils/Utils.jl")
    include("FolsomData/FolsomData.jl")
    include("BegData/BegData.jl")
    include("iJR904/iJR904.jl")

    function __init__()
        _make_dirs()
    end

end
