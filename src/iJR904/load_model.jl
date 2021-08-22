## -------------------------------------------------------------------
function load_model(modelkey::String; uncompress = true)
    file = procdir(iJR904, "base_models.bson")
    models = ldat(file; verbose = false);
    model = models[modelkey]
    return uncompress ? ChU.uncompressed_model(model) : model
end

function load_model(modelkey::String, exp::Int; uncompress = true)
    file = procdir(iJR904, "base_models.bson")
    models = ldat(file; verbose = false);
    model = models[modelkey][exp]
    return uncompress ? ChU.uncompressed_model(model) : model
end