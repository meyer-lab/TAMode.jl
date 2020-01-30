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
    return parse.(Float64, Array(conc)), parse.(Float64, tps), parse.(Float64, Matrix(measVal))
end


@model BLI(tps, conc, bindData) = begin
    Tshift = 822.2
    Kon ~ Normal(6., 0.5)
    Kdis ~ Normal(2., 1.)
    Tshift = tps[1] + 599.9
    theor_save = []
    exp_save = []

    for i in 1:length(conc)
        tBind = tps[tps.<Tshift].-tps[1]
            
        bind_step = R1Calc(conc[i], Kon, Kdis, tBind)
        unbind_step = R2Calc(bind_step[end], Kdis, tps[tps .> Tshift] .- Tshift)
        theor_bind = vcat(bind_step[:], unbind_step)

        if i == 1
            theor_save = theor_bind
            exp_save = bindData[:,i]
        else
            theor_save = vcat(theor_save, theor_bind)
            exp_save = vcat(exp_save, bindData[:,i])
        end
    end

    Rmax = mean(theor_save)/mean(exp_save)
    residuals = exp_save .- (theor_save*Rmax)

    residuals ~ MvNormal(zeros(length(residuals)), ones(length(residuals)) * std(residuals))
end
