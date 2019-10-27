module TAMode

using OrdinaryDiffEq
using StaticArrays

include("reactCode.jl")


function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


function runTAM(tps::Array{Float64,1}, params::Vector)::Array{Float64,2}
    @assert all(params .>= 0.0)
    @assert all(tps .>= 0.0)

    u0 = zeros(55)

    prob = ODEProblem(TAM_reacti, u0, (0.0, maximum(tps)), params)

    sol = solve(prob, Rodas4P(); isoutofdomain=domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, Nspecies))
    end

    return solut
end


end # module
