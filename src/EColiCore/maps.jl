# Here I include some maps between the experimental
# data ids and the model ids

## -------------------------------------------------------------------
# maps between Heerden2013 https://doi.org/10.1186/1475-2859-12-80 
# data and the model


function load_mets_map() 
    Fd_mets_map = Dict()
    Fd_mets_map["GLC"] = "glc__D_e"
    Fd_mets_map["AC"] = "ac_e"
    Fd_mets_map["LAC"] = "lac__D_e"
    Fd_mets_map["PYR"] = "pyr_e"
    Fd_mets_map["SUCC"] = "succ_e"
    Fd_mets_map["FORM"] = "for_e"
    Fd_mets_map["O2"] = "o2_e"
    Fd_mets_map["CO2"] = "co2_e"
    for (k, v) in Fd_mets_map
        Fd_mets_map[v] = k
    end
    return Fd_mets_map
end

# EX_glc_LPAREN_e_RPAREN_
function load_Fd_rxn_map()
    Fd_rxns_map = Dict()
    Fd_rxns_map["D"] = "BIOMASS_Ecoli_core_w_GAM"
    Fd_rxns_map["GLC"] = "EX_glc__D_e"
    Fd_rxns_map["AC"] = "EX_ac_e"
    Fd_rxns_map["LAC"] = "EX_lac__D_e"
    Fd_rxns_map["PYR"] = "EX_pyr_e"
    Fd_rxns_map["SUCC"] = "EX_succ_e"
    Fd_rxns_map["FORM"] = "EX_for_e"
    Fd_rxns_map["O2"] = "EX_o2_e"
    Fd_rxns_map["CO2"] = "EX_co2_e"
    for (k, v) in Fd_rxns_map
        Fd_rxns_map[v] = k
    end
    return Fd_rxns_map
end

# base model exch met map
# A quick way to get exchages from mets and te other way around
function load_exch_met_map()
    !isfile(EXCH_MET_MAP_FILE) && return Dict()
    return load_data(EXCH_MET_MAP_FILE; verbose = false)
end
