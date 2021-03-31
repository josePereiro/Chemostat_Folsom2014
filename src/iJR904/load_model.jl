## -------------------------------------------------------------------
function load_model(modelkey::String; uncompress = true)
    BASE_MODELS = ChU.load_data(BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS[modelkey]
    return uncompress ? ChU.uncompressed_model(model) : model
end

function load_model(modelkey::String, exp::Int; uncompress = true)
    BASE_MODELS = ChU.load_data(BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS[modelkey][exp]
    return uncompress ? ChU.uncompressed_model(model) : model
end