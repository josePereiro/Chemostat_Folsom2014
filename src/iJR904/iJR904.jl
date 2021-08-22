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

    using ProjAssistant
    @gen_sub_proj

    include("const.jl")
    include("load_data.jl")
    include("beg_enz_cost.jl")
    include("load_model.jl")
    include("ME_MODES.jl")
    
    function __init__()
        @create_proj_dirs
    end

end