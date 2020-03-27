const fraction = 0.1 # Fraction of the cell surface with PS


" Setup the parameters for the full TAM receptor model. "
function compParamm(params::Vector{T})::comprates{T} where {T}
    @assert all(params .>= 0.0)

    TAM, hetRs = paramTAMrate(view(params, 4:11))

    out = comprates{T}(TAM, params[1], params[2], params[3], hetRs)

    return detailedBalance(out)
end


function TAMreact(dxxdt_d::Vector, xx_d::Vector, p::comprates, t)
    # Note that TAM_reacti will set ydot to zero first
    # Reaction for the rich phase
    p.rr.gasCur *= p.gasPart
    TAM_reacti(dxxdt_d, xx_d, p.rr, t)

    # Reaction for the poor phase
    p.rr.gasCur /= p.gasPart
    TAM_reacti(dxxdt_d[(Nspecies + 1):end], xx_d[(Nspecies + 1):end], p.rr, t)

    # Deal with transport across boundary
    transIn = p.partIn * p.fraction
    dxxdt_d[1:Nspecies] += (transIn * xx_d[(Nspecies + 1):(2 * Nspecies)] - xx_d[(Nspecies + 1):(2 * Nspecies)] .* boundLig) / p.fraction
    dxxdt_d[(Nspecies + 1):(2 * Nspecies)] -= transIn * (xx_d[(Nspecies + 1):(2 * Nspecies)] - xx_d[1:Nspecies] .* boundLig) / (1 - p.fraction)
end
