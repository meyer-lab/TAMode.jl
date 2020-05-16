const compSize = 100
const boundary = 10 # Nodes with PS
const compDX = compSize * Float64.(collect(1:(compSize - 2)))
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
        @views for ii = 1:27
            duu = view(du, ii:27:length(du))
            uu = u[ii:27:end]

            if boundLigC[ii] == 0
                @. duu[2:(end-1)] = (-4.0*uu[2:(end-1)] + (2.0-1.0/compDX)*uu[1:(end-2)] + (2.0+1.0/compDX)*uu[3:end])/2/dRdRMaxRMaxR
            else
                @. duu[2:(boundary-1)] = (-4.0*uu[2:(boundary-1)] + (2.0-1.0/compDX[1:(boundary-2)])*uu[1:(boundary-2)] + (2.0+1.0/compDX[1:(boundary-2)])*uu[3:boundary])/2/dRdRMaxRMaxR
                duu[boundary] = -4.0*(uu[boundary] - uu[boundary - 1])/dRdRMaxRMaxR

                duu[boundary + 1] = 4.0*(uu[boundary + 2] - uu[boundary + 1])/dRdRMaxRMaxR
                @. duu[(boundary + 2):(end-1)] = (-4.0*uu[(boundary + 2):(end-1)] + (2.0-1.0/compDX[(boundary + 1):end])*uu[(boundary + 1):(end-2)] + (2.0+1.0/compDX[(boundary + 1):end])*uu[(boundary + 3):end])/2/dRdRMaxRMaxR

                # Diffusion into compartment
                duu[boundary + 1] -= uu[boundary + 1] / dRdRMaxRMaxR
                duu[boundary] += uu[boundary + 1] / dRdRMaxRMaxR
            end

            duu[1] = 4.0*(uu[2] - uu[1])/dRdRMaxRMaxR
            duu[end] = -4.0*(uu[end] - uu[end - 1])/dRdRMaxRMaxR
        end

        rmul!(du, r.diff)
    else
        fill!(du, 0.0)
    end

    if reaction
        @views for ii = 0:(Int(length(u) / 27) - 1)
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
