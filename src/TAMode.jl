module TAMode

using OrdinaryDiffEq
using StaticArrays
using SteadyStateDiffEq
using LinearAlgebra
using LabelledArrays
using ForwardDiff

numType = Union{Float64, ForwardDiff.Dual}

include("reactCode.jl")
include("compModel.jl")
include("bothLigands.jl")


function domainDef(u, p, t)
    return any(x -> x < 0.0, u)
end


function getAutocrine(params)
    # TODO: Replace with steady-state
    uNot = convert(Vector{eltype(params)}, zeros(55))
    tps = convert(eltype(params), 10000000.0)

    probInit = ODEProblem(TAM_reacti, uNot, tps, params)
    solInit = solve(probInit, AutoTsit5(TRBDF2()); isoutofdomain=domainDef)

    return solInit(tps)
end


function runTAMinit(tps::Vector, params::Vector, solInit::Vector)
    ansType = promote_type(eltype(tps), eltype(params), eltype(solInit))
    tps = convert(Vector{ansType}, tps)
    solInit = convert(Vector{ansType}, solInit)

    prob = ODEProblem(TAM_reacti, solInit, maximum(tps), params)

    sol = solve(prob, AutoTsit5(TRBDF2()); isoutofdomain=domainDef)
    solut = sol(tps).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end



function runTAM(tps::Array{Float64,1}, params, gasStim::Float64)::Array{Float64,2}
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
