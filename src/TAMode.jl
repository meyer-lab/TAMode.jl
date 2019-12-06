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


function getAutocrine(params::Union{Vector{T}, TAMode.Rates{T}})::Vector{T} where {T}
    probInit = SteadyStateProblem(TAM_reacti, zeros(T, 55), params)

    solInit = solve(probInit, DynamicSS(Rosenbrock23(autodiff=(T == Float64))); isoutofdomain=domainDef)
  
    return solInit.u
end


function runTAMinit(tps::Vector{Float64}, params::Union{Vector{T}, TAMode.Rates{T}}, solInit::Vector) where {T}
    solInit = convert(Vector{T}, solInit)
    prob = ODEProblem(TAM_reacti, solInit, maximum(tps), params)

    sol = solve(prob, Rosenbrock23(autodiff=(T == Float64)); isoutofdomain=domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end



function runTAM(tps::Vector{Float64}, params, gasStim::Float64)
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params)

    if params isa Rates
        params.gasCur = gasStim
    else
        params[7] = gasStim
    end

    return runTAMinit(tps, params, solInit)
end


end # module
