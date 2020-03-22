module TAMode

using OrdinaryDiffEq
using StaticArrays
using SteadyStateDiffEq
using LinearAlgebra
import LabelledArrays
using Turing
import CSV
import Statistics

include("reactCode.jl")
include("bothLigands.jl")
include("compModel.jl")
include("BLI.jl")


const solTol = 1.0e-7

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end

const options = Dict([:reltol => solTol, :abstol => solTol, :isoutofdomain => domainDef])


function getAutocrine(params::Union{Rates{T}, comprates{T}, Lsrates{T}})::Vector{T} where {T <: Real}
    if params isa Rates
        N = 55
    elseif params isa comprates
        N = 110
    elseif params isa Lsrates
        N = 30
    end

    probInit = SteadyStateProblem(TAMreact, zeros(T, N), params)

    sol = AutoTsit5(Rodas5(autodiff = (T == Float64)))
    return solve(probInit, DynamicSS(sol); options...).u
end


function runTAMinit(tps::AbstractVector{Float64}, params::Union{Rates{T}, comprates{T}, Lsrates{T}}, solInit::Vector)::Matrix{T} where {T <: Real}
    @assert all(x -> x >= 0.0, tps)

    solInit = convert(Vector{T}, solInit)
    prob = ODEProblem(TAMreact, solInit, maximum(tps), params)

    sol = AutoTsit5(Rodas5(autodiff = (T == Float64)))
    solut = solve(prob, sol; saveat = tps, options...).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end


function runTAM(tps::AbstractVector{Float64}, params::Union{Rates{T}, Vector{T}}, gasStim::Float64)::Matrix{T} where {T <: Real}
    if params isa Vector
        params = param(params)
    end

    solInit = getAutocrine(params)

    params.gasCur = gasStim

    return runTAMinit(tps, params, solInit)
end


function calcStim(tps::AbstractVector{Float64}, params, gasStim::Float64)
    if params isa Vector
        params = compParamm(params)
    end

    solInit = getAutocrine(params)

    params.rr.gasCur = gasStim

    return runTAMinit(tps, params, solInit)
end


function calcStimPtdser(tps::AbstractVector{Float64}, params)
    solInit = getAutocrine(params)

    return runTAMinit(tps, params, [solInit solInit])
end


function runTAMLS(tps::AbstractVector{Float64}, pIn, ligStim::Tuple{Real, Real})
    params = Lsparam(pIn)

    solInit = getAutocrine(params)
    params.curL = ligStim

    return runTAMinit(tps, params, solInit)
end


include("fitting.jl")

end # module
