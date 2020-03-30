const fraction = 0.1 # Fraction of the cell surface with PS


" Setup the parameters for the full TAM receptor model. "
function compParamm(params::Vector{T})::comprates{T} where {T}
    @assert all(params .>= 0.0)

    TAM, hetRs = paramTAMrate(view(params, 5:12))

    out = comprates{T}(TAM, params[1], params[2], params[3], params[4], hetRs)

    return detailedBalance(out)
end


function TAMreact(du::Vector, u::Vector, r::comprates, t)
    fill!(du, 0.0)
    cache = Vector{promote_type(eltype(du), typeof(r.gasCur))}(undef, 10)

    compartmentReact(u, du, r.gasCur, nothing, r, cache)
end


function TAMreactComp(du::Vector, u::Vector, r::comprates, t)
    fill!(du, 0.0)
    cache = Vector{promote_type(eltype(du), typeof(r.xFwd))}(undef, 10)

    sizze = Int(length(u) / 27)
    boundary = Int(floor(sizze / 10))

    for ii = 0:boundary
        compartmentReact(view(u, (27 * ii + 1):(27 * ii + 27)), view(du, (27 * ii + 1):(27 * ii + 27)), r.gasCur * r.gasPart, nothing, r, cache)
    end

    for ii = boundary:(sizze + 1)
        compartmentReact(view(u, (27 * ii + 1):(27 * ii + 27)), view(du, (27 * ii + 1):(27 * ii + 27)), r.gasCur, nothing, r, cache)
    end

    dx = collect(1:sizze)
    Δ = CenteredDifference(2, 2, dx, sizze, coeff_func = r.diff)
    Δin = CenteredDifference(2, 2, dx[1:boundary], boundary, coeff_func = r.diff)
    Δout = CenteredDifference(2, 2, dx[(boundary + 1):sizze], sizze - boundary, coeff_func = r.diff)
    bc = Neumann0BC(dx, approximation_order)
    bcin = Neumann0BC(dx[1:boundary], approximation_order)
    bcout = Neumann0BC(dx[(boundary + 1):sizze], approximation_order)

    # TODO: Handle boundary flux
    for ii = 1:27
        duu = @view du[ii:27:end]
        uu = u[ii:27:end]

        if boundLigC[ii] == 0
            duu += (Δ * bc) * uu
        else
            duu[1:boundary] += (Δin * bcin) * uu[1:boundary]
            duu[(boundary + 1):end] += (Δout * bcout) * uu[(boundary + 1):end]

            fluxx = (uu[boundary + 1] - uu[boundary]) / sizze / sizze * r.diff
            du[boundary] += fluxx
            du[boundary + 1] -= fluxx
        end
    end
end
