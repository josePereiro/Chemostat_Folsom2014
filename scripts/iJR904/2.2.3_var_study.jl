let
    method = ME_MAX_POL

    var_avs = Dict()
    var_stds = Dict()
    for exp in EXPS
        datfile = INDEX[method, :DFILE, exp]
        dat = deserialize(datfile)
        model = dat[:model]
        idxs = map(Fd.msd_mets) do Fd_ider
            ChU.rxnindex(model, Fd_rxns_map[Fd_ider])
        end
        
        epouts = dat[:epouts]
        exp_beta = maximum(keys(epouts))
        epout = epouts[exp_beta]
        
        vas = ChU.va(epout)[idxs]
        var_avs[exp] = mean(vas)
        var_stds[exp] = 0.0
    end

    # sGLC
    ids = [:D, :X, :uAC, :uGLC, :uLAC, :uPYR]
    for id in ids
        p = scatter(;xlabel = string(id), ylabel = "avs var")
        xs = [Fd.val(id, exp) for exp in EXPS]
        ys = [var_avs[exp] for exp in EXPS]
        errs = [var_stds[exp] for exp in EXPS]
        sids = sortperm(xs)
        xs, ys, errs = xs[sids], ys[sids], errs[sids]
        scatter!(p, xs, ys; yerr = errs, label = "", color = :black, m = 8)
        plot!(p, xs, ys; label = "", color = :black, lw = 3, alpha = 0.5)
        mysavefig(p, string("av_var_vs_", id))
    end
end