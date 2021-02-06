module Chemostat_Folsom2014

    import BSON
    import DrWatson
    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP

    include("Utils/Utils.jl")

    function __init__()
        _make_dirs()
    end

end
