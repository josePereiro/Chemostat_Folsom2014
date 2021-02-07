import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    # -------------------------------------------------------------------
    import Chemostat_Folsom2014
    const ChF = Chemostat_Folsom2014

    const iJR = ChF.iJR904
    const Fd = ChF.FolsomData # experimental data
    const Bd = ChF.BegData    # cost data

    #  ----------------------------------------------------------------------------
    # run add "https://github.com/josePereiro/Chemostat" in the 
    # julia Pkg REPL for installing the package
    import Chemostat
    import Chemostat.LP.MathProgBase
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    import JuMP, GLPK
    import JuMP.MathOptInterface
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
    
    import FileIO
    using Plots
    import GR
    GR.inline("png")
end

## -----------------------------------------------------------------------------------------------
LPDAT = ChU.load_data(iJR.LP_DAT_FILE)
const FBA_BOUNDED = :FBA_BOUNDEDs
const FBA_OPEN = :FBA_OPEN
const YIELD = :YIELD

## -----------------------------------------------------------------------------------------------
fileid = "3.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), iJR.MODEL_FIGURES_DIR; params...)
    @info "Plotting" fname
end
myminmax(a::Vector) = (minimum(a), maximum(a))  
FLX_IDERS = ["GLC", "PYR", "SUCC", "LAC", "FORM", "AC", "O2", "CO2"]

EXPS = 1:4 # experiments that have both concentration and flx data

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = let
    iders = [FLX_IDERS; "D"]
    colors = Plots.distinguishable_colors(length(iders))
    Dict(ider => color for (ider, color) in zip(iders, colors))
end

method_colors = Dict(
    FBA_OPEN => :red,
    FBA_BOUNDED => :orange,
    YIELD => :blue,
)

## -----------------------------------------------------------------------------------------------
# yield correlation
let
    p = plot(;title = "Yield correlation", xlabel = "exp", ylabel = "model")
    m, M = Inf, -Inf
    for exp in EXPS
        try
            model = LPDAT[YIELD, :model, exp]
            yout = LPDAT[YIELD, :yout, exp]
            model_yield = LPDAT[YIELD, :yield, exp]
            
            objidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
            exp_growth = Fd.val("D", exp)
            model_growth = ChU.av(model, yout, objidx)

            diff = abs(model_growth - exp_growth)/exp_growth
            diff > 0.05 && continue # unfeasible

            exp_yield = abs(Fd.val("D", exp) / Fd.uval("GLC", exp))
            scatter!(p, [exp_yield], [model_yield]; ms = 8,
                    color = :blue, alpha = 0.6, label = ""
            )
            m = minimum([m, exp_yield, model_yield])
            M = maximum([M, exp_yield, model_yield])
        catch err; @warn("Fail", err) end
    end
    plot!(p, [m,M], [m,M]; ls = :dash, color = :black, label = "")
    pname = "yield_corr"
    mysavefig(p, pname)
end

## -----------------------------------------------------------------------------------------------
# yield vs stuff
let
    ps = Plots.Plot[]
    for id in [:D, :cGLC, :xi, :uGLC]
        p = plot(;title = "yield vs $(id)", xlabel = "exp $id", ylabel = "yield")
        for exp in EXPS
            try
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]
                model_yield = LPDAT[YIELD, :yield, exp]
                # status, yflxs, model_yield, d, model = DAT[D]
                exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
                biomass_idx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
                
                Fd_val = Fd.val(id, exp)
                exp_yield = abs(Fd.val("D", exp) / Fd.uval("GLC", exp))
                scatter!(p, [Fd_val], [model_yield]; ms = 8,
                        color = :blue, alpha = 0.6, label = ""
                )
                scatter!(p, [Fd_val], [exp_yield]; ms = 8,
                        color = :red, alpha = 0.6, label = ""
                )
            catch err; @warn("Fail", err) end
        end
        push!(ps, p)
    end
    pname = "yield_vs_stuff"
    mysavefig(ps, pname)
end

## -----------------------------------------------------------------------------------------------
iJR.load_exch_met_map()["for_e"]
## -----------------------------------------------------------------------------------------------
# correlations
FLX_IDERS_MAP = Dict(
    "GLC" => "EX_glc_LPAREN_e_RPAREN__REV",
    "AC" => "EX_ac_LPAREN_e_RPAREN_",
    "LAC" => "EX_lac_D_LPAREN_e_RPAREN_",
    "PYR" => "EX_pyr_LPAREN_e_RPAREN_",
    "SUCC" => "EX_succ_LPAREN_e_RPAREN_",
    "FORM" => "EX_for_LPAREN_e_RPAREN_",
    "O2" => "EX_o2_LPAREN_e_RPAREN__REV",
    "CO2" => "EX_co2_LPAREN_e_RPAREN_",
    "D" => iJR.BIOMASS_IDER,
)
    
## -----------------------------------------------------------------------------------------------
# flx correlations
let
    yield_p = plot(title = "yield tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    open_fba_p = plot(title = "open fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    bounded_fba_p = plot(title = "bounded fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    margin, m, M = -Inf, Inf, -Inf
    for (Fd_ider, model_ider) in FLX_IDERS_MAP
        Fd_fun = Fd_ider == "D" ? Fd.val : Fd.uval
        Fd_errf = Fd_ider == "D" ? Fd.err : Fd.uerr
        for exp in EXPS

                color = ider_colors[Fd_ider]
                Fd_flx = abs(Fd_fun(Fd_ider, exp)) 
                Fd_err = abs(Fd_errf(Fd_ider, exp)) 
                
                # yield
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]

                ymax_flx = ChU.av(model, yout, model_ider)
                diffsign = sign(Fd_flx) * sign(ymax_flx)
                Fd_vals = abs(Fd_flx) * diffsign
                ep_vals = abs(ymax_flx) * diffsign

                scatter!(yield_p, [Fd_flx], [ymax_flx]; ms = 8,
                    xerr = [Fd_err],
                    color, alpha = 0.6, label = ""
                )

                # bounded fba
                for (fba_type, p) in [(FBA_BOUNDED, bounded_fba_p) , 
                                    (FBA_OPEN, open_fba_p)]

                    model = LPDAT[fba_type, :model, exp]
                    fbaout = LPDAT[fba_type, :fbaout, exp]
                    
                    fba_flx = ChU.av(model, fbaout, model_ider)
                    scatter!(p, [Fd_flx], [fba_flx]; ms = 8,
                        xerr = [Fd_err],
                        color, alpha = 0.6, label = ""
                    )
                    m = minimum([m, Fd_flx, ymax_flx, fba_flx])
                    M = maximum([M, Fd_flx, ymax_flx, fba_flx])
                end

        end
    end
    margin = abs(M - m) * 0.1
    ps = [yield_p, bounded_fba_p, open_fba_p]
    for p in ps
        plot!(p, [m - margin, M + margin], [m - margin, M + margin]; 
            ls = :dash, color = :black, label = "")
    end
    
    pname = "flx_tot_corr"
    layout = (1, 3)
    mysavefig(ps, pname; layout)
end

## -------------------------------------------------------------------
# join flx correlations
let
    figdir = iJR.MODEL_FIGURES_DIR
    for (lp_p, ep_p, join_name) in [
        ("3.1_flx_tot_corr.png", "2.1_flx_tot_corr.png", "flx_join_corr.png")
    ] 
        lp_img = FileIO.load(joinpath(figdir, lp_p))
        ep_img = FileIO.load(joinpath(figdir, ep_p))
        join_p = UJL.make_grid([lp_img, ep_img])
        fname = joinpath(figdir, join_name)
        FileIO.save(fname, join_p)
        @info "Plotting" fname
    end
end

## -------------------------------------------------------------------
# # # leyends
# # # TODO fix this...
# # let
# #     for (title, colors) in [
# #             ("exp", exp_colors), 
# #             ("iders", ider_colors),
# #             ("method", method_colors)
# #         ]
# #     p = plot(; framestyle = :none)
# #         scolors = sort(collect(colors); by = (p) -> string(first(p)))
# #         for (id, color) in scolors
# #             scatter!(p, [0], [0];
# #                 thickness_scaling = 1,
# #                 color, ms = 8, label = string(id),
# #                 legendfontsize=10, 
# #                 # size = [300, 900],
# #                 # legend = :left
# #             )
# #         end
# #         mysavefig(p, "$(title)_color_legend")
# #     end
# # end