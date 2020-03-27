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

    for ii in 0:boundary
        compartmentReact(view(u, (27*ii + 1):(27*ii + 27)), view(du, (27*ii + 1):(27*ii + 27)), r.gasCur * r.gasPart, nothing, r, cache)
    end

    for ii in boundary:(sizze + 1)
        compartmentReact(view(u, (27*ii + 1):(27*ii + 27)), view(du, (27*ii + 1):(27*ii + 27)), r.gasCur, nothing, r, cache)
    end
    
    dx = collect(1:sizze)
    Δ = CenteredDifference(2, 2, dx, sizze, coeff_func=r.diff)

    # TODO: Handle boundary conditions
    for ii in 1:27
        if boundLigC[ii] == 0
            du[ii:27:end] += Δ * u[ii:27:end]
        else
            du[ii:27:end] += Δ * u[ii:27:end]
        end
    end
end
