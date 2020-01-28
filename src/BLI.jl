
""" Calculation for binding step. """
function R1Calc(conc::Real, Kon::Real, Kdis::Real, tps)
    return conc / (Kdis / Kon + conc) * (1 .- (1 / (exp.((Kon * conc + Kdis) * tps))))
end


""" Calculation for unbinding step. """
function R2Calc(R::Real, Kdis::Real, tps)
    return R * exp.(-Kdis * tps)
end
