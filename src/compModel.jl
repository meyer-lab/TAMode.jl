
mutable struct comprates{T}
    rr::Rates{T}
    fraction::T #Fraction of cell surface covered with PtdSer
    partIn::T #Partitioning rate into PtdSer regions
    gasPart::T #Partitioning of ligand into PtdSer region
end


function compParamm(compIn::Vector)
    return comprates{eltype(compIn)}(param(compIn[4:end]), compIn[1], compIn[2], compIn[3])
end


function getAutocrineComp(params::Union{Vector{T}, comprates{T}})::Vector{T} where {T}
    probInit = SteadyStateProblem(TAMreactComp, zeros(T, 110), params)

    sol = Rodas5(autodiff = (T == Float64))
    return solve(probInit, DynamicSS(sol); options...).u
end


function TAMreactComp(dxxdt_d, xx_d, p, t)
    if p isa Vector
        p = compParamm(p)
    end

    # Note that TAM_reacti will set ydot to zero first
    # Reaction for the rich phase
    p.rr.gasCur *= p.gasPart
    TAM_reacti(dxxdt_d, xx_d, p.rr, t) 

    # Reaction for the poor phase
    p.rr.gasCur /= p.gasPart
    TAM_reacti(dxxdt_d[Nspecies + 1:end], xx_d[Nspecies + 1:end], p.rr, t)

    # Deal with transport across boundary
    transIn = p.partIn * p.fraction
    dxxdt_d[1:Nspecies] += (transIn*xx_d[Nspecies + 1 : 2*Nspecies] - xx_d[Nspecies + 1 : 2*Nspecies] .* boundLig) / p.fraction
    dxxdt_d[Nspecies + 1 : 2*Nspecies] -= transIn*(xx_d[Nspecies + 1 : 2*Nspecies] - xx_d[1:Nspecies] .* boundLig) / (1 - p.fraction)
end


function calcStim(tps::Array{Float64,1}, params, gasStim::Float64)
    @assert all(tps .>= 0.0)

    solInit = getAutocrineComp(params)

    if params isa Rates
        params.gasCur = gasStim
    else
        params[7] = gasStim
    end

    prob = ODEProblem(TAMreactComp, solInit, maximum(tps), params)
    solut = solve(prob, Rodas5(); saveat = tps, options...).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end


function calcStimPtdser(tps::Array{Float64,1}, params)
    @assert all(tps .>= 0.0)

    solInit = getAutocrine(params, TAM_reacti, 55)
    solInit = [solInit solInit]

    prob = ODEProblem(TAMreactComp, solInit, maximum(tps), params)
    solut = solve(prob, Rodas5(); saveat = tps, options...).u

    if length(tps) > 1
        solut = vcat(transpose.(solut)...)
    else
        solut = reshape(solut[1], (1, length(solInit)))
    end

    return solut
end

