import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    ## -------------------------------------------------------------------
    import CSV
    import MAT
    using DataFrames
    using Serialization

    ## -------------------------------------------------------------------
    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014
    const Fd  = ChF.FolsomData;
    const ECC = ChF.EColiCore

    ## -------------------------------------------------------------------
    import Chemostat
    const ChU = Chemostat.Utils
    const ChSS = Chemostat.SteadyState
    const ChLP = Chemostat.LP
end

## -------------------------------------------------------------------
# Tools
function partial_test(model, exp = 4; title  = "PARTIAL TEST")
    fbaout = ChLP.fba(model, ECC.BIOMASS_IDER);
    ChU.tagprintln_inmw(title, 
        "\nsize:             ", size(model),
        "\nobj_ider:         " , ECC.BIOMASS_IDER,
        "\nfba obj_val:      ", ChU.av(model, fbaout, ECC.BIOMASS_IDER),
        "\nexp obj_val:      ", Fd.val("D", exp),
    )
end

## -------------------------------------------------------------------
# BASE MODEL
# Load Mat file
println("Original .mat model")
src_file = ECC.MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["e_coli_core"]
model = ChU.MetNet(mat_model; reshape = true);

## -------------------------------------------------------------------
# the maximal experimental growth rate in FOlsom2014 is ~0.4 1/h
# The raw model present a growth rate bigger than that, so it is ok
# to use it directly as base model
partial_test(model)

## -------------------------------------------------------------------
# Set bounds
# The abs maximum bounds will be set to ECC.MAX_ABS_BOUND
ChU.tagprintln_inmw("CLAMP BOUNDS", 
    "\nabs max bound: ", ECC.MAX_ABS_BOUND
)
ChU.clampfields!(model, [:lb, :ub]; 
    abs_max = ECC.MAX_ABS_BOUND, zeroth = 0.0)
partial_test(model)

## -------------------------------------------------------------------
# Exchanges
exchs = ChU.exchanges(model)
ChU.tagprintln_inmw("SETTING EXCHANGES") 
# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with xi = 1.0
# see Cossios paper (see README)
ξ = minimum(Fd.val(:xi))
# println("Minimum medium: ", ECC.base_intake_info)
foreach(exchs) do idx
    ChU.ub!(model, idx, ECC.MAX_ABS_BOUND) # Opening all outakes
    ChU.lb!(model, idx, 0.0) # Closing all intakes
end

# see Cossios paper (see README) for details in the Chemostat bound constraint
intake_info = ECC.load_base_intake_info()
ChSS.apply_bound!(model, ξ, intake_info);
partial_test(model)

## -------------------------------------------------------------------
# Eliminate similars
ChU.tagprintln_inmw("ELIMINATE SIMILAR RXNS") 
ChU.summary.([model], ["SUCDi", "FRD7"])
ChU.bounds!(model, "SUCDi", 0.0, 0.0)
ChU.bounds!(model, "FRD7", -ECC.MAX_ABS_BOUND, ECC.MAX_ABS_BOUND)
ChU.summary.([model], ["SUCDi", "FRD7"])

## -------------------------------------------------------------------
# Exch_met_map

# A quick way to get exchages from mets and te other way around
exch_met_map = Dict()
exchis = findall((id) -> any(startswith.(id, ["EX_", "DM_"])), model.rxns)
for exch_ in model.rxns[exchis]
    mets_ = model.mets[ChU.rxn_mets(model, exch_)]
    length(mets_) != 1 && continue
    exch_met_map[exch_] = mets_[1]
    exch_met_map[mets_[1]] = exch_
end;

# saving
ChU.save_data(ECC.EXCH_MET_MAP_FILE, exch_met_map)

## -------------------------------------------------------------------
MODEL_DAT = Dict()
MODEL_DAT[:load_model] = deepcopy(model) |> ChU.compressed_model;

## -------------------------------------------------------------------
# FVA PREPROCESSING
ChU.tagprintln_inmw("DOING FVA PREPROCESS", "\n")

for (exp, xi) in enumerate(Fd.val(:xi))

    DAT = get!(MODEL_DAT, :fva_model, Dict())
    model0 = deepcopy(model)
    
    exp_intake_info = ECC.load_base_intake_info(exp)
    ChSS.apply_bound!(model0, xi, exp_intake_info);

    fva_model = ChLP.fva_preprocess(model0, 
        batchlen = 1,
        check_obj = ECC.BIOMASS_IDER,
        verbose = true
    )
    partial_test(fva_model, exp)
    DAT[exp] = fva_model |> ChU.compressed_model
    println()
end

## -------------------------------------------------------------------
# MAX POL
let
    # This model is bounded by the maximum rates found for EColi.
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # Extract max exchages from FIG 3 to form the maximum polytope

    max_model = deepcopy(model)
    
    # Biomass
    # 2.2 1/ h
    ChU.bounds!(max_model, ECC.BIOMASS_IDER, 0.0, 2.2)
    
    Fd_rxn_map = ECC.load_Fd_rxn_map()
    # 40 mmol / gDW h
    ChU.bounds!(max_model, Fd_rxn_map["GLC"], -40.0, 0.0)
    # Reduce box to improve convergency
    # ChU.bounds!(max_model, Fd_rxn_map["GLC"], -10.0, 0.0)
    # 45 mmol/ gDW
    ChU.bounds!(max_model, Fd_rxn_map["AC"], 0.0, 40.0)
    # 55 mmol/ gDW h
    ChU.bounds!(max_model, Fd_rxn_map["FORM"], 0.0, 55.0)
    # 20 mmol/ gDW h
    ChU.bounds!(max_model, Fd_rxn_map["O2"], -20.0, 0.0)
    
    # fva
    max_model = ChLP.fva_preprocess(max_model, 
        check_obj = ECC.BIOMASS_IDER,
        verbose = true
    );
    
    partial_test(max_model)

    ## -------------------------------------------------------------------
    # saving
    MODEL_DAT[:max_model] = ChU.compressed_model(max_model)
end;


 ## -------------------------------------------------------------------
# Saving model
ChU.save_data(ECC.BASE_MODELS_FILE, MODEL_DAT);