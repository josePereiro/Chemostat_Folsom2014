module FolsomData

    include("data.jl")

    function __init__()
        _populate_bundle()
    end

end