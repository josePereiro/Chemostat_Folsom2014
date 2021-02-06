import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Kayser2005")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
    # Run add https://github.com/josePereiro/Chemostat_Kayser2005.jl in the Julia Pkg REPL to install the
    # package, then you must activate the package enviroment (see README)
    import Chemostat_Kayser2005
    const ChK = Chemostat_Kayser2005

    const iJR = ChK.iJR904
    const Kd = ChK.KayserData # experimental data
    const Bd = ChK.BegData    # cost data

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
CONC_IDERS = ["GLC", "AC", "NH4"]
FLX_IDERS = ["GLC", "CO2", "O2", "AC", "NH4"]

EXPS = 1:13 # experiments that have both concentration and flx data

exp_colors = let
    colors = Plots.distinguishable_colors(length(EXPS))
    Dict(exp => color for (exp, color) in zip(EXPS, colors))
end

ider_colors = Dict(
    "GLC" => :red, "CO2" => :yellow,
    "O2" => :blue, "AC" => :orange, 
    "NH4" => :green, "D" => :black,
)

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
            
            objidx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
            exp_growth = Kd.val("D", exp)
            model_growth = ChU.av(model, yout, objidx)

            diff = abs(model_growth - exp_growth)/exp_growth
            diff > 0.05 && continue # unfeasible

            exp_yield = abs(Kd.val("D", exp) / Kd.uval("GLC", exp))
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
    for id in [:D, :cGLC, :xi, :uGLC, :sGLC]
        p = plot(;title = "yield vs $(id)", xlabel = "exp $id", ylabel = "yield")
        for exp in EXPS
            try
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]
                model_yield = LPDAT[YIELD, :yield, exp]
                # status, yflxs, model_yield, d, model = DAT[D]
                exglc_idx = ChU.rxnindex(model, "EX_glc_LPAREN_e_RPAREN__REV")
                biomass_idx = ChU.rxnindex(model, iJR.KAYSER_BIOMASS_IDER)
                
                Kd_val = Kd.val(id, exp)
                exp_yield = abs(Kd.val("D", exp) / Kd.uval("GLC", exp))
                scatter!(p, [Kd_val], [model_yield]; ms = 8,
                        color = :blue, alpha = 0.6, label = ""
                )
                scatter!(p, [Kd_val], [exp_yield]; ms = 8,
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
# correlations
FLX_IDERS_MAP = Dict(
    "GLC" => "EX_glc_LPAREN_e_RPAREN__REV",
    "CO2" => "EX_co2_LPAREN_e_RPAREN_",
    "O2" => "EX_o2_LPAREN_e_RPAREN__REV",
    "AC" => "EX_ac_LPAREN_e_RPAREN_",
    "NH4" => "EX_nh4_LPAREN_e_RPAREN__REV",
    "D" => iJR.KAYSER_BIOMASS_IDER
)
    
SENSE = Dict(
    "GLC" => -1,
    "CO2" => 1,
    "O2" => -1,
    "AC" => 1,
    "NH4" => -1,
    "D" => 1
)
## -----------------------------------------------------------------------------------------------
# flx correlations
let
    yield_p = plot(title = "yield tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    open_fba_p = plot(title = "open fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    bounded_fba_p = plot(title = "bounded fba tot corrs"; xlabel = "exp flx", ylabel = "model flx")
    margin, m, M = -Inf, Inf, -Inf
    for (Kd_ider, model_ider) in FLX_IDERS_MAP
        Kd_fun = Kd_ider == "D" ? Kd.val : Kd.uval
        for exp in EXPS

                color = ider_colors[Kd_ider]
                Kd_flx = abs(Kd_fun(Kd_ider, exp)) # every body is possitive here
                
                # yield
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]

                ymax_flx = ChU.av(model, yout, model_ider)
                diffsign = sign(Kd_flx) * sign(ymax_flx)
                Kd_vals = abs(Kd_flx) * diffsign
                ep_vals = abs(ymax_flx) * diffsign

                scatter!(yield_p, [Kd_flx], [ymax_flx]; ms = 8,
                    color, alpha = 0.6, label = ""
                )

                # bounded fba
                for (fba_type, p) in [(FBA_BOUNDED, bounded_fba_p) , 
                                    (FBA_OPEN, open_fba_p)]

                    model = LPDAT[fba_type, :model, exp]
                    fbaout = LPDAT[fba_type, :fbaout, exp]
                    
                    fba_flx = ChU.av(model, fbaout, model_ider)
                    scatter!(p, [Kd_flx], [fba_flx]; ms = 8,
                        color, alpha = 0.6, label = ""
                    )
                    m = minimum([m, Kd_flx, ymax_flx, fba_flx])
                    M = maximum([M, Kd_flx, ymax_flx, fba_flx])
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

## -----------------------------------------------------------------------------------------------
# conc correlations
let
    yield_p = plot(title = "yield tot corrs"; xlabel = "exp conc", ylabel = "model conc")
    open_fba_p = plot(title = "open fba tot corrs"; xlabel = "exp conc", ylabel = "model conc")
    bounded_fba_p = plot(title = "bounded fba tot corrs"; xlabel = "exp conc", ylabel = "model conc")
    
    margin, m, M = -Inf, Inf, -Inf
    for (Kd_ider, model_ider) in FLX_IDERS_MAP
        Kd_ider in ["D", "CO2", "O2"] && continue
        for exp in EXPS

                color = ider_colors[Kd_ider]
                Kd_sval = Kd.sval(Kd_ider, exp)
                Kd_cval = Kd.cval(Kd_ider, exp, 0.0)
                exp_xi = Kd.val(:xi, exp)

                # yield
                !haskey(LPDAT, YIELD, :model, exp) && continue
                model = LPDAT[YIELD, :model, exp]
                yout = LPDAT[YIELD, :yout, exp]

                # conc (s = c + u*xi)
                ymax_flx = ChU.av(model, yout, model_ider)
                ymax_sval =  max(Kd_cval + SENSE[Kd_ider] * ymax_flx * exp_xi, 0.0)

                scatter!(yield_p, [Kd_sval], [ymax_sval]; ms = 8,
                    color, alpha = 0.6, label = ""
                )

                # bounded fba
                for (fba_type, p) in [(FBA_BOUNDED, bounded_fba_p) , 
                                    (FBA_OPEN, open_fba_p)]

                    model = LPDAT[fba_type, :model, exp]
                    fbaout = LPDAT[fba_type, :fbaout, exp]
                    
                    # conc (s = c + u*xi)
                    fba_flx = ChU.av(model, fbaout, model_ider)
                    fba_sval = Kd_cval == 0 ? fba_flx * exp_xi :
                        max(Kd_cval - fba_flx * exp_xi, 0.0)
                    scatter!(p, [Kd_sval], [fba_sval]; ms = 8,
                        color, alpha = 0.6, label = ""
                    )
                    m = minimum([m, Kd_sval, ymax_sval, fba_sval])
                    M = maximum([M, Kd_sval, ymax_sval, fba_sval])
                end

        end
    end
    margin = abs(M - m) * 0.1
    ps = [yield_p, bounded_fba_p, open_fba_p]
    for p in ps
        plot!(p, [m - margin, M + margin], [m - margin, M + margin]; 
            ls = :dash, color = :black, label = "")
    end
    
    pname = "conc_tot_corr"
    layout = (1, 3)
    mysavefig(ps, pname; layout)
end

## -------------------------------------------------------------------
# join flx correlations
let
    figdir = iJR.MODEL_FIGURES_DIR
    for (lp_p, ep_p, join_name) in [
        ("3.1_flx_tot_corr.png", "2.1_flx_tot_corr.png", "flx_join_corr.png"),
        ("3.1_conc_tot_corr.png", "2.1_conc_tot_corr.png", "conc_join_corr.png"),
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