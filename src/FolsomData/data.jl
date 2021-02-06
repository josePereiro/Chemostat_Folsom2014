## ----------------------------------------------------------------------------
# Data from the glucose-limited culture reported at
# Folsom, James Patrick, Albert E. Parker, and Ross P. Carlson. 
# “Physiological and Proteomic Analysis of Escherichia Coli Iron-Limited 
# Chemostat Growth.” Journal of Bacteriology 196, no. 15 (August 1, 2014): 
# 2748–61. https://doi.org/10.1128/JB.01606-14.

## ----------------------------------------------------------------------------
# Sumpplementary table 2
TABLE_S2 = Dict()
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
TABLE_S2["qglc"] = Dict(
    "unit" => "Cmmol/gCDW-h",
    "val" => [9.1, 14.8, 25.5, 33.9],
    "err" => [0.68, 0.78, 1.87, 2.79],
)

## ----------------------------------------------------------------------------
#  M9 minima medium, all in mM
MEDIUM = Dict()
MEDIUM["glc"] = 2.2 # mM


## ----------------------------------------------------------------------------
# Bundle, converting
DATA = Dict()
let
    DATA["D"] = deepcopy(TABLE_S2["D"])
    DATA["μ"] = deepcopy(TABLE_S2["μ"])

    DATA["qglc"] = Dict(
        # q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "unit" => "mmol/gCDW-h",
        "val" => TABLE_S2["qglc"]["val"] ./ 6,
        "err" => TABLE_S2["qglc"]["err"] ./ 6,
    )
    
    qglc = TABLE_S2["qglc"]["val"]
    DATA["pyr"] = Dict(
        # Y[Cmmol/Cmmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Ypyr_glc"]["val"] ./ qglc ./ 3,
        "err" => TABLE_S2["Ypyr_glc"]["err"] ./ qglc ./ 3,
    )
    DATA["succ"] = Dict(
        # Y[Cmmol/Cmmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Ysucc_glc"]["val"] ./ qglc ./ 4,
        "err" => TABLE_S2["Ysucc_glc"]["err"] ./ qglc ./ 4,
    )
    DATA["lac"] = Dict(
        # Y[Cmmol/Cmmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Ylac_glc"]["val"] ./ qglc ./ 3,
        "err" => TABLE_S2["Ylac_glc"]["err"] ./ qglc ./ 3,
    )
    DATA["form"] = Dict(
        # Y[Cmmol/Cmmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Yform_glc"]["val"] ./ qglc ./ 1,
        "err" => TABLE_S2["Yform_glc"]["err"] ./ qglc ./ 1,
    )
    DATA["ac"] = Dict(
        # Y[Cmmol/Cmmol] * qglc[Cmmol/gCDW-h] = q[Cmmol/gCDW-h] / #C = q[mmol/gCDW-h]
        "unit" => "mmol/gCDW h",
        "val" => TABLE_S2["Yac_glc"]["val"] ./ qglc ./ 2,
        "err" => TABLE_S2["Yac_glc"]["err"] ./ qglc ./ 2,
    )

    # There negigeble glucose present in the medium at steady state 
    # so we can compute X = c * D/qglc
    DATA["X"] = Dict(
        # c[mmol/L] * D[1/h] / q[mmol/gCDW h] = X[gCDW/L]
        "val" => MEDIUM["glc"] .* DATA["D"]["val"] ./ DATA["qglc"]["val"],
        "err" => zero(DATA["qglc"]["err"]) #TODO: fix this
    )

    # xi
    DATA["xi"] = Dict(
        "val" => DATA["X"]["val"] ./ DATA["D"]["val"],
        "err" => DATA["X"]["err"] #TODO: fix this
    )

end

# ## ------------------------------------------------------------------
# # API
# const EXPS = 1:4
# const msd_mets = ["AC", "GLC", "NH4"]
# const iders_to_plot = ["AC", "GLC", "NH4", "D"]
# val(dataid) = BUNDLE[string(dataid)]["val"]
# val(dataid, exp::Int) = val(dataid)[exp]
# function val(dataid, exp)
#     dat = BUNDLE[string(dataid)]
#     i = findfirst(dat["D"] .== exp)
#     val(dataid)[i]
# end
# val(dataid, exp, dflt) = try; return val(dataid, exp); catch err; (@warn(err); dflt) end

# cval(id, args...) = val("c$id", args...)
# sval(id, args...) = val("s$id", args...)
# uval(id, args...) = val("u$id", args...)

# name(dataid) = BUNDLE[string(dataid)]["name"]
# unit(dataid) = BUNDLE[string(dataid)]["unit"]

