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
    using Serialization
    import UtilsJL
    const UJL = UtilsJL
end

## -----------------------------------------------------------------------------------------------
const FBA_BOUNDED = :FBA_BOUNDEDs
const FBA_OPEN = :FBA_OPEN
const YIELD = :YIELD
const EXPS = 1:13

## -----------------------------------------------------------------------------------------------
function base_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model_dict = BASE_MODELS["fva_models"][exp]
    model = ChU.MetNet(;model_dict...) |> ChU.uncompressed_model
end

## -----------------------------------------------------------------------------------------------

# for method in [HOMO, EXPECTED, BOUNDED]
LPDAT = UJL.DictTree()
let
    objider = iJR.KAYSER_BIOMASS_IDER
    glcider = "EX_glc_LPAREN_e_RPAREN__REV"
    iterator = Kd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator
        try 
            # prepare model
            model = base_model(exp)
            M, N = size(model)
            model = ChU.fix_dims(model)
            ChU.invert_bkwds!(model)
            exglc_idx = ChU.rxnindex(model, glcider)
            objidx = ChU.rxnindex(model, objider)
            exp_growth = Kd.val(:D, exp)
            dgrowth = 0.01
            
            ChU.ub!(model, objider, exp_growth * (1 + dgrowth))

            # fba
            fbaout = ChLP.fba(model, objider)
            fba_growth = ChU.av(model, fbaout, objider)

            # yield max
            model.c .= 0.0; model.c[objidx] = 1.0
            d = zeros(N); d[exglc_idx] = 1.0
            status, yflxs, yield = ChLP.yLP(model, d)
            status != JuMP.MOI.OPTIMAL && @warn(status)
            ymax_growth = yflxs[objidx]
            LPDAT[YIELD, :status, exp] = status
            LPDAT[YIELD, :yield, exp] = yield
            LPDAT[YIELD, :yout, exp] = ChU.FBAout(yflxs, yflxs[objidx], objider, nothing)
            LPDAT[YIELD, :d, exp] = d
            LPDAT[YIELD, :model, exp] = model

            @info("Yield Maximization", 
                exp, status, yield,
                fba_growth, ymax_growth, exp_growth
            ); println()

        catch err
            @warn("Error", err, exp); println()
        end
    end
end

## -------------------------------------------------------------------
# FBA
let
    objider = iJR.KAYSER_BIOMASS_IDER
    costider = iJR.COST_IDER
    biomass_f = 0.01

    iterator = Kd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator

        @info("Doing ", exp); println()

        dat = Dict()

        # OPEN
        model = base_model(exp)
        model = ChU.fix_dims(model)
        ChU.invert_bkwds!(model)
        fbaout = ChLP.fba(model, objider, costider)
        fba_growth = ChU.av(model, fbaout, objider)
        
        # storing
        LPDAT[FBA_OPEN, :model, exp] = model
        LPDAT[FBA_OPEN, :fbaout, exp] = fbaout

        # BOUNDED
        model = base_model(exp)
        model = ChU.fix_dims(model)
        ChU.invert_bkwds!(model)
        objidx = ChU.rxnindex(model, objider)
        exp_growth = Kd.val("D", exp)
        ChU.ub!(model, objider, exp_growth * (1.0 + biomass_f))
        fbaout = ChLP.fba(model, objider, costider)
        fba_growth = ChU.av(model, fbaout, objider)

        # storing
        LPDAT[FBA_BOUNDED, :model, exp] = model
        LPDAT[FBA_BOUNDED, :fbaout, exp] = fbaout
    end
end

## -------------------------------------------------------------------
ChU.save_data(iJR.LP_DAT_FILE, LPDAT)
