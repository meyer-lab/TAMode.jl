module TAMode

using OrdinaryDiffEq
using StaticArrays
using SteadyStateDiffEq
using LinearAlgebra

include("reactCode.jl")
include("compModel.jl")


function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


function getAutocrine(params::Vector, funcc, nZero::Int)
    @assert all(params .>= 0.0)
    
    # TODO: Replace with steady-state
    probInit = ODEProblem(TAM_reacti, zeros(nZero), 10000000.0, params)
    solInit = solve(probInit, AutoTsit5(Rosenbrock23()))
    
    return solInit(10000000.0)
end


function runTAM(tps::Array{Float64,1}, params::Vector, gasStim::Float64)::Array{Float64,2}
    @assert all(params .>= 0.0)
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params, TAM_reacti, 55)

    params[7] = gasStim
    prob = ODEProblem(TAM_reacti, solInit, maximum(tps), params)

    sol = solve(prob, AutoTsit5(Rosenbrock23()); isoutofdomain=domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end


end # module
