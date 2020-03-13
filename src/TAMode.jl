module TAMode

using OrdinaryDiffEq
using StaticArrays
using SteadyStateDiffEq
using LinearAlgebra
using LabelledArrays
using Turing
using CSV
using Statistics

include("reactCode.jl")
include("bothLigands.jl")
include("compModel.jl")
include("BLI.jl")


const solTol = 1.0e-9

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end

const options = Dict([:reltol => solTol, :abstol => solTol, :isoutofdomain => domainDef])


function getAutocrine(params::Union{Vector{T}, Rates{T}, comprates{T}, Lsrates{T}}, func, N::Int)::Vector{T} where {T}
    probInit = SteadyStateProblem(func, zeros(T, N), params)

    sol = Rodas5(autodiff = (T == Float64))
    return solve(probInit, DynamicSS(sol); options...).u
end


function runTAMinit(tps::Vector{Float64}, params::Union{Vector{T}, Rates{T}, comprates{T}, Lsrates{T}}, func, solInit::Vector) where {T}
    solInit = convert(Vector{T}, solInit)
    prob = ODEProblem(func, solInit, maximum(tps), params)

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

    solInit = getAutocrine(params, TAM_reacti, 55)

    if params isa Rates
        params.gasCur = gasStim
    else
        params[7] = gasStim
    end

    return runTAMinit(tps, params, TAM_reacti, solInit)
end


function calcStim(tps::Vector{Float64}, params, gasStim::Float64)
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params, TAMreactComp, 110)

    if params isa comprates
        params.gasCur = gasStim
    else
        params[7] = gasStim
    end

    return runTAMinit(tps, params, TAMreactComp, solInit)
end


function calcStimPtdser(tps::Vector{Float64}, params)
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params, TAM_reacti, 55)

    return runTAMinit(tps, params, TAMreactComp, [solInit solInit])
end


function runTAMLS(tps::Vector{Float64}, pIn, ligStim::Tuple{Real, Real})
    params = Lsparam(pIn)
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params, TAMreactLS, 30)
    params.curL = ligStim

    return runTAMinit(tps, params, TAMreactLS, solInit)
end


include("fitting.jl")

end # module
