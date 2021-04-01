# # DIRS
# const PROJ_NAME = "Folsom2014"
# const FOLSOM_PROCESSED_DATA_DIR = joinpath(PROCESSED_DATA_DIR, PROJ_NAME)
# const FOLSOM_FIGURES_DIR = joinpath(FIGURES_DATA_DIR, PROJ_NAME)

# function _create_dirs()
#     for dir in [FOLSOM_PROCESSED_DATA_DIR, FOLSOM_FIGURES_DIR]
#         try; mkpath(dir); catch err; @warn("Error creating dir", dir, err); end
#     end
# end