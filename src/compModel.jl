
mutable struct comprates{T} #do we want this to be T and what does {T} mean
    rr::Rates #do i need values for all the rates and does this want to be in the reactCode file or a new one
    fraction::T #Fraction of cell surface covered with PtdSer
    partIn::T #Partitioning rate into PtdSer regions
    gasPart::T #Partitioning of ligand into PtdSer region
end

function calcStim(tps::Array{Float64,1}, params, gasStim::Float64)::Array{Float64,2}
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


function calcStimPtdser(tps::Array{Float64,1}, params, gasStim::Float64)::Array{Float64,2}
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

