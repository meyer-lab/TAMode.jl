const fraction = 0.1 # Fraction of the cell surface with PS


" Setup the parameters for the full TAM receptor model. "
function compParamm(params::Vector{T})::comprates{T} where {T}
    @assert all(params .>= 0.0)

    TAM, hetRs = paramTAMrate(view(params, 5:12))

    out = comprates{T}(TAM, params[1], params[2], params[3], params[4], hetRs)

    return detailedBalance(out)
end


function getDiffOp(dx, D)
    Δ = CenteredDifference(1, 2, dx, length(dx) - 1, D)
    return Δ * Neumann0BC(dx, 1)
end


function TAMreact(du::Vector, u::Vector, r::comprates, t; reaction = true)
    fill!(du, 0.0)

    sizze = Int(length(u) / 27)
    boundary = Int(floor(sizze / 10))

    if reaction
        for ii = 0:(sizze - 1)
            if ii < boundary
                gass = r.gasCur * r.gasPart
            else
                gass = r.gasCur
            end

            idx = (27 * ii + 1):(27 * ii + 27)
            compartmentReact(view(u, idx), view(du, idx), gass, nothing, r)
        end
    end

    # If we're not yet dealing with the PDE, solve for starting state
    if length(du) == 27
        return nothing
    end

    dx = Float64.(collect(1:(sizze + 1)))
    bc = getDiffOp(dx, r.diff)
    bcin = getDiffOp(dx[1:(boundary + 1)], r.diff)
    bcout = getDiffOp(dx[boundary:sizze], r.diff)

    for ii = 1:27
        duu = @view du[ii:27:end]
        uu = @view u[ii:27:end]

        if boundLigC[ii] == 0
            duu += bc * uu
        else
            duu[1:boundary] += bcin * uu[1:boundary]
            duu[(boundary + 1):end] += bcout * uu[(boundary + 1):end]
            # TODO: Implement one-way flux across boundary
        end
    end
end
