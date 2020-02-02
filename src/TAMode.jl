module TAMode

using OrdinaryDiffEq
using StaticArrays
using SteadyStateDiffEq
using LinearAlgebra
using LabelledArrays

include("reactCode.jl")
include("bothLigands.jl")


const solTol = 1.0e-9

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end

const options = Dict([:reltol => solTol, :abstol => solTol, :isoutofdomain => domainDef])


function getAutocrine(params::Union{Vector{T}, TAMode.Rates{T}})::Vector{T} where {T}
    probInit = SteadyStateProblem(TAM_reacti, zeros(T, 55), params)

    sol = Rodas5(autodiff = (T == Float64))
    return solve(probInit, DynamicSS(sol); options...).u
end


function runTAMinit(tps::Vector{Float64}, params::Union{Vector{T}, TAMode.Rates{T}}, solInit::Vector) where {T}
    solInit = convert(Vector{T}, solInit)
    prob = ODEProblem(TAM_reacti, solInit, maximum(tps), params)

    sol = Rodas5(autodiff = (T == Float64))
    solut = solve(prob, sol; saveat = tps, options...).u

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


include("compModel.jl")


end # module
