module EColiCore

import ..Chemostat_Folsom2014.FolsomData
const Fd = FolsomData
import ..Chemostat_Folsom2014: PROJ_ROOT, DATA_DIR, FIGURES_DATA_DIR, RAW_DATA_DIR, PROCESSED_DATA_DIR
import Chemostat.Utils: load_data

include("const.jl")
include("file_and_dirs.jl")
include("base_intake_info.jl")
include("maps.jl")

function __init__()
    _create_dirs()
end

end