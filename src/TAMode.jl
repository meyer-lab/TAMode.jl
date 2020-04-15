module TAMode

using OrdinaryDiffEq
using StaticArrays
using LinearAlgebra
import LabelledArrays
using Turing
import Statistics
using DiffEqOperators
using ModelingToolkit

include("types.jl")
include("reactCode.jl")
include("bothLigands.jl")
include("compModel.jl")


const solTol = 1.0e-9

function domainDef(u, p, t)
    return any(x -> x < -solTol, u)
end


function getAutocrine(params::Union{Rates{T}, comprates{T}, Lsrates{T}})::Vector{T} where {T <: Real}
    if params isa Rates
        N = 55
    elseif params isa comprates
        N = 27
    elseif params isa Lsrates
        N = 30
    end

    return vec(runTAMinit([1.0e6], params, zeros(T, N)))
end


function runTAMinit(tps::AbstractVector{Float64}, params::Union{Rates{T}, comprates{T}, Lsrates{T}}, solInit::Vector{T})::Matrix{T} where {T <: Real}
    @assert all(x -> x >= 0.0, tps)

    prob = ODEProblem(TAMreact, solInit, maximum(tps), params)

    sol = AutoTsit5(Rodas5(autodiff = (T == Float64)))
    solut = solve(prob, sol; saveat = tps, reltol = solTol, isoutofdomain = domainDef).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end


function compTAM(tps::AbstractVector{Float64}, params::Union{comprates{T}, Vector{T}}) where {T <: Real}
    if params isa Vector
        params = compParamm(params)
    end

    solInit = getAutocrine(params)
    u0 = repeat(solInit; outer = [100])

    return runTAMinit(tps, params, u0)
end


function runTAM(tps::AbstractVector{Float64}, params::Union{Rates{T}, Lsrates{T}, Vector{T}}, ligStim)::Matrix{T} where {T <: Real}
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
    elseif params isa Lsrates
        params.curL = ligStim
    else
        @assert false
    end

    return runTAMinit(tps, params, solInit)
end


include("fitting.jl")

end # module
