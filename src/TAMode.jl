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


const solTol = 1.0e-5

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


function runTAM(tps::AbstractVector{Float64}, params::Union{Rates{T}, comprates{T}, Lsrates{T}, Vector{T}}, ligStim)::Matrix{T} where {T <: Real}
    if params isa Vector
        if length(params) == 15
            params = param(params)
        elseif length(params) == 9
            params = Lsparam(params)
        else
            @assert false
        end
    end

    solInit = getAutocrine(params)

    if params isa Rates
        params.gasCur = ligStim
    elseif params isa comprates
        params.rr.gasCur = ligStim
    elseif params isa Lsrates
        params.curL = ligStim
    else
        @assert false
    end

    return runTAMinit(tps, params, solInit)
end


include("fitting.jl")

end # module
