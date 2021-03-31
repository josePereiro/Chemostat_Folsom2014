import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import MAT

    import SparseArrays

    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    const iJR = ChF.iJR904
    const Fd = ChF.FolsomData # experimental data
    const Bd = ChF.BegData    # cost data

    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP

end

## ------------------------------------------------------------------
# LOAD RAW MODEL
src_file = iJR.MODEL_RAW_MAT_FILE
mat_model = MAT.matread(src_file)["model"]
model = ChU.MetNet(mat_model; reshape=true)
ChU.tagprintln_inmw("MAT MODEL LOADED", 
    "\nfile:             ", relpath(src_file), 
    "\nfile size:        ", filesize(src_file), " bytes", 
    "\nmodel size:       ", size(model),
    "\nChU.nzabs_range:      ", ChU.nzabs_range(model.S),
)
ChF.test_fba(model, iJR.BIOMASS_IDER; summary = false)

## -------------------------------------------------------------------
# Set bounds
# The abs maximum bounds will be set to 100
ChU.tagprintln_inmw("CLAMP BOUNDS", 
    "\nabs max bound: ", iJR.ABS_MAX_BOUND
)
foreach(model.rxns) do ider
        ChU.isfixxed(model, ider) && return # fixxed reaction are untouched

        old_ub = ChU.ub(model, ider)
        new_ub = old_ub == 0.0 ? 0.0 : iJR.ABS_MAX_BOUND
        ChU.ub!(model, ider, new_ub)

        old_lb = ChU.lb(model, ider)
        new_lb = old_lb == 0.0 ? 0.0 : -iJR.ABS_MAX_BOUND
        ChU.lb!(model, ider, new_lb)
end

## -------------------------------------------------------------------
# CLOSING EXCHANGES
exchs = ChU.exchanges(model)
ChU.tagprintln_inmw("CLOSE EXCANGES", 
    "\nChU.exchanges: ", exchs |> length
)
# Close, for now, all ChU.exchanges for avoiding it to be in revs
# The reversible reactions will be splited for modeling cost
# Exchanges have not associated cost, so, we do not split them
foreach(exchs) do idx
    ChU.ub!(model, idx, 0.0) # Closing all outtakes
    ChU.lb!(model, idx, 0.0) # Closing all intakes
end

## -------------------------------------------------------------------
# # EXCHANGE METABOLITE MAP
# exch_met_map = Dict()
# for rxni in exchs
#     rxn = model.rxns[rxni]
#     metis = ChU.rxn_mets(model, rxni)
#     mets = model.mets[metis]
#     length(mets) != 1 && continue 
#     met = mets |> first
#     exch_met_map[rxn] = met
#     exch_met_map[met] = rxn
# end
# ChU.save_data(iJR.EXCH_MET_MAP_FILE, exch_met_map)

## -------------------------------------------------------------------
# ENZYMATIC COST INFO
# The cost will be introduced as a reaction, we follow the same cost models as 
# Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
# A new balance equations is then added:
#        Σ(rᵢ*costᵢ) + tot_cost = 0
#    Because the cost coefficients (costᵢ) < 0 (it resamble a reactant), the system must allocate 
#    the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
#    are usually bounded [0.0, 1.0]
cost_info = Dict()
fwd_ider(rxn) = string(rxn, ChU.FWD_SUFFIX);
bkwd_ider(rxn) = string(rxn, ChU.BKWD_SUFFIX);
for rxn in model.rxns
    # The ChU.exchanges, the atpm and the biomass are synthetic reactions, so, 
    # they have should not have an associated enzimatic cost 
    any(startswith.(rxn, ["EX_", "DM_"])) && continue
    rxn == iJR.BIOMASS_IDER && continue
    rxn == iJR.BIOMASS_IDER && continue
    rxn == iJR.ATPM_IDER && continue
        
    # Only the internal, non reversible reactions have an associated cost
    # We will split the rev reactions, so we save the cost for both versions (fwd, bkwd)
    if ChU.isrev(model, rxn)
        cost_info[fwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
        cost_info[bkwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
    else
        cost_info[rxn] = -iJR.beg_enz_cost(rxn)
    end
end

## -------------------------------------------------------------------
# SPLITING REVS
ChU.tagprintln_inmw("SPLITING REVS", 
    "\nfwd_suffix:      ", ChU.FWD_SUFFIX,
    "\nbkwd_suffix:     ", ChU.BKWD_SUFFIX,
)
model = ChU.split_revs(model;
    get_fwd_ider = fwd_ider,
    get_bkwd_ider = bkwd_ider,
);

## -------------------------------------------------------------------
# ADDING COST REACCION
cost_met_id = "cost"
cost_exch_id = iJR.COST_IDER
ChU.tagprintln_inmw("ADDING COST", 
    "\ncosts to add: ", cost_info |> length,
    "\nmin abs coe:  ", cost_info |> values .|> abs |> minimum,
    "\nmax abs coe:  ", cost_info |> values .|> abs |> maximum,
    "\ncost met id:  ", cost_met_id,
    "\ncost exch id: ", cost_exch_id
)

M, N = size(model)
cost_met = ChU.Met(cost_met_id, S = collect(values(cost_info)), rxns = collect(keys(cost_info)), b = 0.0)
model = ChU.expanded_model(model, M + 1, N + 1)
ChU.set_met!(model, ChU.findempty(model, :mets), cost_met)
cost_exch = ChU.Rxn(cost_exch_id, S = [1.0], mets = [cost_met_id], lb = -iJR.ABS_MAX_BOUND, ub = 0.0, c = 0.0)
ChU.set_rxn!(model, ChU.findempty(model, :rxns), cost_exch);

## -------------------------------------------------------------------
# SET BASE EXCHANGE
ChU.tagprintln_inmw("SETTING EXCHANGES") 
# To control the intakes just the metabolites defined in the 
# base_intake_info (The minimum medium) will be opened.
# The base model will be constraint as in a cultivation with 
# experimental minimum xi
# see Cossios paper (see README)

foreach(exchs) do idx
    ChU.ub!(model, idx, iJR.ABS_MAX_BOUND) # Opening all outakes
    ChU.lb!(model, idx, 0.0) # Closing all intakes
end

# see Cossios paper (see README) for details in the Chemostat bound constraint
xi = minimum(Fd.val(:xi))
intake_info = iJR.load_base_intake_info()
ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

# tot_cost is the exchange that controls the bounds of the 
# enzimatic cost contraint, we bound it to [0, 1.0]
ChU.lb!(model, cost_exch_id, 0.0);
ChU.ub!(model, cost_exch_id, 1.0);

## -------------------------------------------------------------------
model = ChU.fix_dims(model)
ChF.test_fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER)

## -------------------------------------------------------------------
# FVA PREPROCESSING
const BASE_MODELS = isfile(iJR.BASE_MODELS_FILE) ? 
    ChU.load_data(iJR.BASE_MODELS_FILE) : 
    Dict{Any, Any}("load_model" => ChU.compressed_model(model))

## -------------------------------------------------------------------
let
    for (exp, D) in Fd.val(:D) |> enumerate

        DAT = get!(BASE_MODELS, "fva_models", Dict())
        ChU.tagprintln_inmw("DOING FVA", 
            "\nexp:             ", exp,
            "\nD:               ", D,
            "\ncProgress:       ", length(DAT),
            "\n"
        )
        haskey(DAT, exp) && continue # cached

        ## -------------------------------------------------------------------
        # prepare model
        model0 = deepcopy(model)
        M, N = size(model0)
        exp_xi = Fd.val(:xi, exp)
        intake_info = iJR.intake_info(exp)
        ChSS.apply_bound!(model0, exp_xi, intake_info; 
            emptyfirst = true)

        ChF.test_fba(exp, model0, iJR.BIOMASS_IDER, iJR.COST_IDER)
        fva_model = ChLP.fva_preprocess(model0, 
            check_obj = iJR.BIOMASS_IDER,
            verbose = true
        );
        ChF.test_fba(exp, fva_model, iJR.BIOMASS_IDER, iJR.COST_IDER)

        # storing
        DAT[exp] = ChU.compressed_model(fva_model)

        ## -------------------------------------------------------------------
        # caching
        ChU.save_data(iJR.BASE_MODELS_FILE, BASE_MODELS);
        GC.gc()
    end
end

## -------------------------------------------------------------------
# MAX MODEL
let
    # This model is bounded by the maximum rates found for EColi.
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # Extract max exchages from FIG 3 to form the maximum polytope

    ChU.tagprintln_inmw("DOING MAX MODEL", "\n")

    max_model = deepcopy(model)
    
    # Biomass
    # 2.2 1/ h
    ChU.bounds!(max_model, iJR.BIOMASS_IDER, 0.0, 2.2)
    
    Fd_rxns_map = iJR.load_Fd_rxns_map() 
    # 40 mmol / gDW h
    ChU.bounds!(max_model, Fd_rxns_map["GLC"], -40.0, 0.0)
    # 45 mmol/ gDW
    ChU.bounds!(max_model, Fd_rxns_map["AC"], 0.0, 40.0)
    # 55 mmol/ gDW h
    ChU.bounds!(max_model, Fd_rxns_map["FORM"], 0.0, 55.0)
    # 20 mmol/ gDW h
    ChU.bounds!(max_model, Fd_rxns_map["O2"], -20.0, 0.0)
    
    # fva
    max_model = ChLP.fva_preprocess(max_model, 
        check_obj = iJR.BIOMASS_IDER,
        verbose = true
    );

    ## -------------------------------------------------------------------
    test_model = deepcopy(max_model)
    for exp in 1:4
        D = Fd.val(:D, exp)
        cgD_X = Fd.cval(:GLC, exp) * Fd.val(:D, exp) / Fd.val(:X, exp)
        ChU.lb!(test_model, iJR.GLC_EX_IDER, -cgD_X)
        fbaout = ChLP.fba(test_model, iJR.BIOMASS_IDER, iJR.COST_IDER)
        biom = ChU.av(test_model, fbaout, iJR.BIOMASS_IDER)
        cost = ChU.av(test_model, fbaout, iJR.COST_IDER)
        @info("Test", exp, cgD_X, D, biom, cost); println()
    end

    ## -------------------------------------------------------------------
    # saving
    BASE_MODELS["max_model"] = ChU.compressed_model(max_model)
    ChU.save_data(iJR.BASE_MODELS_FILE, BASE_MODELS);
end;