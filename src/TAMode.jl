module TAMode

using OrdinaryDiffEq
using StaticArrays
using SteadyStateDiffEq
using LinearAlgebra
using LabelledArrays

include("reactCode.jl")
include("compModel.jl")
include("bothLigands.jl")


function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


function getAutocrine(params, funcc, nZero::Int)
    # TODO: Replace with steady-state
    probInit = ODEProblem(TAM_reacti, zeros(nZero), 10000000.0, params)

    solInit = solve(probInit, AutoTsit5(TRBDF2(linsolve=LinSolveGMRES())); isoutofdomain=domainDef)

    return solInit(10000000.0)
end


function runTAM(tps::Array{Float64,1}, params, gasStim::Float64)::Array{Float64,2}
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params, TAM_reacti, 55)

    if params isa Rates
        params.gasCur = gasStim
    else
        params[7] = gasStim
    end

    prob = ODEProblem(TAM_reacti, solInit, maximum(tps), params)

    sol = solve(prob, AutoTsit5(TRBDF2(linsolve=LinSolveGMRES())); isoutofdomain=domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end


end # module
