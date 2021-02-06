# I map some model ids to the enzimatic cost data reported in 
# Beg2007, so, I will set the cost of this reactions using this data
# The rest of the model reactions will get the average enzimatic cost also 
# reported by Beg.

function load_beg_enz_cost_dict()
    beg_enz_cost_dict = Dict()
    beg_rxns_map = load_beg_rxns_map()
    beg_enz_data = Bd.load_enz_data()
    for (i, beg_rxn) in enumerate(beg_enz_data.Enzyme)
        for model_rnx in beg_rxns_map[beg_rxn]
            if model_rnx != ""
                beg_enz_cost_dict[model_rnx] = 
                    beg_enz_data[i, Symbol("ai/xi (h g / mmol)")]
            end
        end
    end
    beg_enz_cost_dict
end

function beg_enz_cost(model_rxn)
    beg_enz_cost_dict = load_beg_enz_cost_dict()
    return get(beg_enz_cost_dict, model_rxn, Bd.AVE_COST)
end