module FolsomData

    import ..Chemostat_Folsom2014: 
        PROCESSED_DATA_DIR, FIGURES_DATA_DIR
    
    include("data.jl")
    include("dir_and_files.jl")

    function __init__()
        _create_dirs()
        _populate_bundle()
    end

end