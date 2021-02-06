const PROJ_ROOT = dirname(dirname(@__DIR__))

const DATA_DIR = joinpath(PROJ_ROOT, "data")

const RAW_DATA_DIR = joinpath(DATA_DIR, "raw")
const PROCESSED_DATA_DIR = joinpath(DATA_DIR, "processed")
const FIGURES_DATA_DIR = joinpath(DATA_DIR, "figures")

function _make_dirs()
    for dir in [DATA_DIR, RAW_DATA_DIR, 
                PROCESSED_DATA_DIR, FIGURES_DATA_DIR]
        if !isdir(dir)
            mkpath(dir)
            println("created $dir")
        end
    end
end