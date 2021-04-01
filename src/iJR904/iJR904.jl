module iJR904

    import BSON
    import ..Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    import ..BegData
    const Bd = BegData
    import ..FolsomData
    const Fd = FolsomData
    import Chemostat
    const ChU = Chemostat.Utils

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)

    include("const.jl")
    include("dirs_and_files.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("load_model.jl")
    function __init__()
        _create_dirs()
        UJL.create_proj_dirs(@__MODULE__)
    end

end