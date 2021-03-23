# The intakes bounds of the network are determined by the 
# medium concentration in the Chemostat model (see Cossios paper)
# This is a base medium for modeling
const MAX_CONC = 99999.0
const ABS_MAX_BOUND = 1000.0

function load_base_intake_info()
    return Dict(
        "EX_glc__D_e" => Dict("c"=> maximum(Fd.val(:cGLC)), "lb"=> -ABS_MAX_BOUND),
        
        "EX_nh4_e" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_o2_e"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_pi_e"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        # "SO4" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND), # Not found in model
    )
end

function load_base_intake_info(exp)
    return Dict(
        "EX_glc__D_e" => Dict("c"=> Fd.val(:cGLC, exp), "lb"=> -ABS_MAX_BOUND),
        "EX_nh4_e" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_o2_e"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        "EX_pi_e"  => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND),
        # "SO4" => Dict("c"=> MAX_CONC, "lb"=> -ABS_MAX_BOUND), # Not found in model
    )
end