function test_fba(model, obj_ider, cost_ider = nothing; summary = false)
    fbaout = isnothing(cost_ider) ? ChLP.fba(model, obj_ider) : ChLP.fba(model, obj_ider, cost_ider)
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", obj_ider,
        "\nfba obj_val:      ", ChU.av(model, fbaout, obj_ider),
        "\nmax exp obj_val:  ", maximum(FolsomData.val("D")),
        "\ncost_ider:        ", isnothing(cost_ider) ? "Not included" : cost_ider,
        "\nfba cost_val:     ", isnothing(cost_ider) ? "Not included" : ChU.av(model, fbaout, cost_ider),
    )
    summary && ChU.summary(model, fbaout)
    nothing
end

function test_fba(exp::Int, model, obj_ider, cost_ider = nothing; summary = false)
    fbaout = isnothing(cost_ider) ? ChLP.fba(model, obj_ider) : ChLP.fba(model, obj_ider, cost_ider)
    ChU.tagprintln_inmw("FBA SOLUTION", 
        "\nobj_ider:         ", obj_ider,
        "\nfba obj_val:      ", ChU.av(model, fbaout, obj_ider),
        "\nexp obj_val:      ", FolsomData.val("D", exp),
        "\ncost_ider:        ", isnothing(cost_ider) ? "Not included" : cost_ider,
        "\nfba cost_val:     ", isnothing(cost_ider) ? "Not included" : ChU.av(model, fbaout, cost_ider),
    )
    summary && ChU.summary(model, fbaout)
    nothing
end