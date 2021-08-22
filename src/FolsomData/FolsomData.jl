module FolsomData

    using ProjAssistant
    @gen_sub_proj
    
    include("data.jl")

    function __init__()
        _populate_bundle()
        @create_proj_dirs
    end

end