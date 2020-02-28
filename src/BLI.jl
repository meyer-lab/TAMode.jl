gas6_T1 = "../data/T1-010820.csv"
gas6_T2 = "../data/T2-010820.csv"
gas6_TFL = "../data/TFL-010820.csv"


""" Calculation for binding step. """
function R1Calc(conc::Real, Kon::Real, Kdis::Real, tps)
    return conc / (Kdis / Kon + conc) * (1 .- (1 ./ (exp.((Kon .* conc .+ Kdis) * tps))))
end


""" Calculation for unbinding step. """
function R2Calc(R::Real, Kdis::Real, tps)
    return R * exp.(-Kdis * tps)
end


function importData(cond)
    df = CSV.read(cond)
    conc = df[1, 2:end]
    tps = df[5:end, 1]
    measVal = df[5:end, 2:end]
    return parse.(Float64, Array(conc)), parse.(Float64, tps), parse.(Float64, Matrix(measVal))
end


function bindingCalc(tps::Vector, Kon::Real, Kdis::Real, Rmax::Real)
    Tshift = tps[1] + 599.9
    tBind = tps[tps .< Tshift] .- tps[1]
    tUnbind = tps[tps .> Tshift] .- Tshift

    bind_step = R1Calc(conc[i], Kon, Kdis, tBind)
    unbind_step = R2Calc(bind_step[end], Kdis, tUnbind)
    theor_bind = vcat(bind_step[:], unbind_step) * Rmax

    return theor_bind
end


@model BLI(tps, conc, bindData) = begin
    Kon ~ LogNormal(6.0, 0.5)
    Kdis ~ LogNormal(1.0, 1.0)
    Rmax ~ LogNormal(-1.0, 0.1)

    resid_save = []

    for i = 1:(length(conc) - 1)
        theor_bind = bindingCalc(tps, Kon, Kdis, Rmax)

        if i == 1
            resid_save = bindData[:, i] .- theor_bind
        else
            resid_save = vcat(resid_save, bindData[:, i] .- theor_bind)
        end
    end

    resid_save ~ MvNormal(zeros(length(resid_save)), ones(length(resid_save)) * std(resid_save))
end
