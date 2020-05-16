const fraction = 0.1 # Fraction of the cell surface with PS
const compSize = 100
const boundary = 10
const compDX = compSize * Float64.(collect(1:(compSize + 1)))
const dRdRMaxRMaxR = 1.0/compSize/compSize


" Setup the parameters for the full TAM receptor model. "
function compParamm(params::Vector{T})::comprates{T} where {T}
    @assert all(params .>= 0.0)

    TAM, hetRs = paramTAMrate(view(params, 5:12))

    out = comprates{T}(TAM, params[1], params[2], params[3], params[4], hetRs)

    return detailedBalance(out)
end


function TAMreact(du::Vector, u::Vector, r::comprates, t; reaction = true)
    # If we're dealing with the PDE form
    if length(du) > 100
        for ii = 1:27
            duu = view(du, ii:27:length(du))
            uu = u[ii:27:end]

            if boundLigC[ii] == 0
                mul!(duu, bc, uu)
            else
                mul!(view(duu, 1:boundary), bcin, uu[1:boundary])
                mul!(view(duu, (boundary + 1):length(duu)), bcout, uu[(boundary + 1):end])
                # TODO: Implement one-way flux across boundary
            end
        end

        rmul!(du, r.diff)
    else
        fill!(du, 0.0)
    end

    if reaction
        for ii = 0:(Int(length(u) / 27) - 1)
            if ii < boundary
                gass = r.gasCur * r.gasPart
            else
                gass = r.gasCur
            end

            idx = (27 * ii + 1):(27 * ii + 27)
            compartmentReact(view(u, idx), view(du, idx), gass, nothing, r)
        end
    end
end
