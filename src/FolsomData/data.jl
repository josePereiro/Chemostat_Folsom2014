## ----------------------------------------------------------------------------
# BUNDLE from the glucose-limited culture reported at
# Folsom, James Patrick, Albert E. Parker, and Ross P. Carlson. 
# “Physiological and Proteomic Analysis of Escherichia Coli Iron-Limited 
# Chemostat Growth.” Journal of Bacteriology 196, no. 15 (August 1, 2014): 
# 2748–61. https://doi.org/10.1128/JB.01606-14.

## ----------------------------------------------------------------------------
# Sumpplementary table 2
TABLE_S2 = Dict()
function _populate_table_s2()
    empty!(TABLE_S2)
    TABLE_S2["D"] = Dict(
        "unit" => "1/h",
        "val" => [0.1, 0.2, 0.3, 0.4],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    TABLE_S2["μ"] = deepcopy(TABLE_S2["D"])
    TABLE_S2["YX_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.426, 0.525, 0.458, 0.459],
        "err" => [0.032, 0.028, 0.034, 0.038],
    )
    TABLE_S2["Ypyr_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.0, 0.0, 0.0, 0.0],
        "err" => [0.0, 0.0, 0.0, 0.0],
    )
    TABLE_S2["Ysucc_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.004, 0.003, 0.0, 0.0],
        "err" => [0.001, 0.001, 0.0, 0.0],
    )
    TABLE_S2["Ylac_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.013, 0.003, 0.026, 0.023],
        "err" => [0.003, 0.001, 0.002, 0.003],
    )
    TABLE_S2["Yform_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.041, 0.005, 0.001, 0.002],
        "err" => [0.015, 0.001, 0.000, 0.002],
    )
    TABLE_S2["Yac_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.000, 0.000, 0.001, 0.004],
        "err" => [0.000, 0.000, 0.001, 0.002],
    )
    TABLE_S2["Yco2_glc"] = Dict(
        "unit" => "Cmmol/Cmmol glc",
        "val" => [0.515, 0.465, 0.514, 0.512],
        "err" => [0.051, 0.031, 0.037, 0.047],
    )

    TABLE_S2["uglc"] = Dict(
        "unit" => "Cmmol/gCDW-h",
        "val" => [9.1, 14.8, 25.5, 33.9],
        "err" => [0.68, 0.78, 1.87, 2.79],
    )
    TABLE_S2["uO2"] = Dict(
        "unit" => "Cmmol/gCDW-h",
        "val" => [4.74, 6.64, 12.68, 16.82],
        "err" => [0.826, 0.786, 1.85, 2.93],
    )
    TABLE_S2
end

## ----------------------------------------------------------------------------
#  M9 minima medium, all in mM
MEDIUM = Dict()
function _populate_medium()
    empty!(MEDIUM)
    MEDIUM["glc"] = 2.2 # mM
    MEDIUM
end

## ----------------------------------------------------------------------------
# Bundle, converting
BUNDLE = Dict()
function _populate_bundle()

    _populate_table_s2()
    _populate_medium()
    empty!(BUNDLE)

    # flxs
    BUNDLE["D"] = deepcopy(TABLE_S2["D"])
    BUNDLE["μ"] = deepcopy(TABLE_S2["μ"])
    BUNDLE["uGLC"] = Dict(
        # q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "glucose uptake rate",
        "unit" => "mmol/gCDW h",
        "val" => -TABLE_S2["uglc"]["val"] ./ 6,
        "err" => TABLE_S2["uglc"]["err"] ./ 6,
    )
    BUNDLE["uO2"] = Dict(
        # q[mmol/gCDW-h]
        "name" => "oxigen uptake rate",
        "unit" => "mmol/gCDW h",
        "val" => -TABLE_S2["uO2"]["val"],
        "err" => TABLE_S2["uO2"]["err"]
    )
    
    # yields
    # qglc[Cmmol/gCDW-h]
    qglc = TABLE_S2["uglc"]["val"]

    BUNDLE["uPYR"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "pyruvate production rate",
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Ypyr_glc"]["val"] .* qglc ./ 3,
        "err" => TABLE_S2["Ypyr_glc"]["err"] .* qglc ./ 3,
    )
    BUNDLE["uSUCC"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "succinate production rate",
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Ysucc_glc"]["val"] .* qglc ./ 4,
        "err" => TABLE_S2["Ysucc_glc"]["err"] .* qglc ./ 4,
    )
    BUNDLE["uLAC"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "lactate production rate",
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Ylac_glc"]["val"] .* qglc ./ 3,
        "err" => TABLE_S2["Ylac_glc"]["err"] .* qglc ./ 3,
    )
    BUNDLE["uFORM"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "formate production rate",
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Yform_glc"]["val"] .* qglc ./ 1,
        "err" => TABLE_S2["Yform_glc"]["err"] .* qglc ./ 1,
    )
    BUNDLE["uAC"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "acetate production rate",
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Yac_glc"]["val"] .* qglc ./ 2,
        "err" => TABLE_S2["Yac_glc"]["err"] .* qglc ./ 2,
    )
    BUNDLE["uCO2"] = Dict(
        # Y[Cmol/Cmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "name" => "acetate production rate",
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Yco2_glc"]["val"] .* qglc ./ 1,
        "err" => TABLE_S2["Yco2_glc"]["err"] .* qglc ./ 1,
    )

    # There negigeble glucose present in the medium at steady state 
    # so we can compute X = c * D/qglc
    BUNDLE["X"] = Dict(
        # c[mmol/L] * D[1/h] / q[mmol/gCDW h] = X[gCDW/L]
        "name" => "cell concentration",
        "unit" => "gCDW/ L",
        "val" => abs.(MEDIUM["glc"] .* BUNDLE["D"]["val"] ./ BUNDLE["uGLC"]["val"]),
        "err" => zero(BUNDLE["uGLC"]["err"]) #TODO: fix this
    )

    # xi
    BUNDLE["xi"] = Dict(
        "name" => "dilution specific cell concentration",
        "unit" => "gCDW/ L h",
        "val" => BUNDLE["X"]["val"] ./ BUNDLE["D"]["val"],
        "err" => BUNDLE["X"]["err"] #TODO: fix this
    )

    # medium
    BUNDLE["cGLC"] = Dict(
        "name" => "glucose feed concentration",
        "unit" => "gCDW/ L h",
        "val" => fill(MEDIUM["glc"], 4),
        "err" => zeros(4)
    )

    BUNDLE
end

## ------------------------------------------------------------------
# API
const EXPS = 1:4
const msd_mets = ["GLC", "PYR", "SUCC", "LAC", "FORM", "AC"]

_get_val(id, dk) = BUNDLE[string(id)][dk]
_get_val(id, dk, exp::Int) = BUNDLE[string(id)][dk][exp]
function _get_val(id, dk, D::Float64)
    exp = findfirst(BUNDLE["D"][dk] .== D)
    isnothing(exp) && error("No experiment with D = $D")
    _get_val(id, dk, exp)
end
_get_val(id, dk, ref, dflt) =
    try; _get_val(id, dk, ref); catch err; dflt end

for fun in [:val, :err]
    dk = string(fun)
    @eval $fun(id) = _get_val(id, $dk)
    @eval $fun(id, ref) = _get_val(id, $dk, ref)
    @eval $fun(id, ref, dflt) = _get_val(id, $dk, ref, dflt)

    for p in [:u, :c]
        pstr = string(p)
        pfun = Symbol(p, fun)
        @eval $pfun(id, args...) = $fun(string($pstr, id), args...)
    end
end

name(id) = BUNDLE[string(id)]["name"]
unit(id) = BUNDLE[string(id)]["unit"]

ciD_X(id) = [cval(id, exp, 0.0) * val(:D, exp) / val(:X, exp) for exp in EXPS]
ciD_X(id, exp) = cval(id, exp, 0.0) * val(:D, exp) / val(:X, exp)