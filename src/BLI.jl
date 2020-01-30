gas6_T1 = "../data/T1-010820.csv"
gas6_T2 = "../data/T2-010820.csv"
gas6_TFL = "../data/TFL-010820.csv"


""" Calculation for binding step. """
function R1Calc(conc::Real, Kon::Real, Kdis::Real, tps)
    return conc / (Kdis / Kon + conc) * (1 .- (1 / (exp.((Kon * conc + Kdis) * tps))))
end


""" Calculation for unbinding step. """
function R2Calc(R::Real, Kdis::Real, tps)
    return R * exp.(-Kdis * tps)
end


function importData(cond)
    df = CSV.read(cond)
    conc = df[1,2:end]
    tps = df[5:end,1]
    measVal = df[5:end,2:end]
    return parse.(Float64,Array(conc)), parse.(Float64, tps), parse.(Float64,Matrix(measVal))
end


@model BLI(tps,conc,bindData) = begin
    Tshift = 822.2
    Kon ~ Normal(6., 0.5)
    Kdis ~ Normal(2., 1.)
    Tshift = tps[1] + 599.9
    theor_save = 0
    exp_save = 0
    for i in 1:length(conc)
        tBind = tps[tps.<Tshift].-tps[1]
        L0 = conc[i]
        exp_bind = bindData[:,i]
            
        bind_step = R1Calc(L0, Kon, Kdis, tBind)
        unbind_step = R2Calc(bind_step[end], Kdis, tps[tps .> Tshift] .- Tshift)
        theor_bind = vcat(bind_step[:], unbind_step)
        if i == 1
            theor_save = theor_bind
            exp_save = exp_bind
        else
            theor_save = vcat(theor_save, theor_bind)
            exp_save = vcat(exp_save, exp_bind)
        end
    end

    Rmax = mean(theor_save)/mean(exp_save)

    residuals = exp_save .- (theor_save*Rmax)
    N = length(residuals)
    #resMean = mean(residuals)
    #resStd = stdm(residuals, resMean)
    #residual ~ Normal(resMean, abs(resStd))
    residuals ~ MvNormal(zeros(N), 2*ones(N))
end
