module FolsomData

    import ..Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    import UtilsJL
    const UJL = UtilsJL
    UJL.gen_sub_proj(@__MODULE__)
    
    
    include("data.jl")
    include("dir_and_files.jl")

    function __init__()
        # _create_dirs()
        _populate_bundle()
        UJL.create_proj_dirs(@__MODULE__)
    end

end