import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Folsom2014")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    #  ----------------------------------------------------------------------------
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
    const ChLP = Ch.LP

    using Serialization
    import UtilsJL
    const UJL = UtilsJL
end

## -----------------------------------------------------------------------------------------------
const FBA_Z_FIX_MIN_COST    = :FBA_Z_FIX_MIN_COST
const FBA_MAX_BIOM_MIN_COST = :FBA_MAX_BIOM_MIN_COST
const FBA_Z_FIX_MIN_VG_COST = :FBA_Z_FIX_MIN_VG_COST
const FBA_Z_VG_FIX_MIN_COST = :FBA_Z_VG_FIX_MIN_COST
const EXPS = 1:4

## -----------------------------------------------------------------------------------------------
function load_model(exp)
    BASE_MODELS = ChU.load_data(iJR.BASE_MODELS_FILE; verbose = false);
    model = BASE_MODELS["fva_models"][exp] |> ChU.uncompressed_model
end

## -----------------------------------------------------------------------------------------------
# Data container
LPDAT = UJL.DictTree()

## -------------------------------------------------------------------
# FBA
let
    objider = iJR.BIOMASS_IDER
    costider = iJR.COST_IDER
    exglcider = iJR.EX_GLC_IDER
    max_sense = -1.0
    min_sense = 1.0

    iterator = Fd.val(:D) |> enumerate |> collect 
    for (exp, D) in iterator

        @info("Doing ", exp); println()

        # FBA_MAX_BIOM_MIN_COST
        let
            model = load_model(exp)
            fbaout = ChLP.fba(model, objider, costider)
            
            LPDAT[FBA_MAX_BIOM_MIN_COST, :model, exp] = model
            LPDAT[FBA_MAX_BIOM_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_COST
        let
            model = load_model(exp)
            exp_growth = Fd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout = ChLP.fba(model, objider, costider)

            LPDAT[FBA_Z_FIX_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MIN_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_FIX_MIN_VG_COST
        let
            model = load_model(exp)
            exp_growth = Fd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            fbaout1 = ChLP.fba(model, exglcider; sense = max_sense)
            exglc = ChU.av(model, fbaout1, exglcider)
            ChU.bounds!(model, exglcider, exglc, exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_FIX_MIN_VG_COST, :model, exp] = model
            LPDAT[FBA_Z_FIX_MIN_VG_COST, :fbaout, exp] = fbaout
        end

        # FBA_Z_VG_FIX_MIN_COST
        let
            model = load_model(exp)
            exp_growth = Fd.val("D", exp)
            ChU.bounds!(model, objider, exp_growth, exp_growth)
            exp_exglc = Fd.uval("GLC", exp)
            ChU.bounds!(model, exglcider, exp_exglc, exp_exglc)
            fbaout = ChLP.fba(model, costider; sense = min_sense)

            LPDAT[FBA_Z_VG_FIX_MIN_COST, :model, exp] = model
            LPDAT[FBA_Z_VG_FIX_MIN_COST, :fbaout, exp] = fbaout
        end
    end
end

## -------------------------------------------------------------------
ChU.save_data(iJR.LP_DAT_FILE, LPDAT)
