Nspecies = 55
partSpecies::Vector{Bool, Nspecies}

mutable struct comprates{T} #do we want this to be T and what does {T} mean
    rr::Rates #do i need values for all the rates and does this want to be in the reactCode file or a new one
    fraction::T #Fraction of cell surface covered with PtdSer
    partIn::T #Partitioning rate into PtdSer regions
    gasPart::T #Partitioning of ligand into PtdSer region
end

function compParamm(compIn::Vector)
    ppparams = comprates{eltype(compIn)}
    ppparams.rr = param(compIn[4:end])
    ppparams.fraction = compIn[1]
    ppparams.partIn = compIn[2]
    ppparams.gasPart = compIn[3]
    
    return ppparams
end

function TAM_react(xx_d, dxxdt_d) #xx_d and dxxdt_d are both pointers, would be just pass an array?
    dxxdt_d = zeros(Float64, Nspecies*2+1)
    p.rr.gasCur *= p.gasPart
    
    TAM_reacti(xx_d,dxxdt_d, p.rr)
     p.rr.gasCur /= p.gasPart        
    
    TAM_reacti(xx_d+Nspecies,dxxdt_d+Nspecies, &p.rr)
        
    transIn = p.partIn*p.fraction
    
    for spec in 1:Nspecies
        if partSpecies[spec]
            dxxdt_d[spec] += transIn*xx_d[Nspecies + spec] / p.fraction
            dxxdt_d[spec + Nspecies] -= transIn*xx_d[Nspecies + spec] / (1 - p.fraction)
        else
            dxxdt_d[spec] += transIn*(xx_d[Nspecies + spec] - xx_d[spec]) / p.fraction
            dxxdt_d[spec + Nspecies] -= transIn*(xx_d[Nspecies + spec] - xx_d[spec]) / (1 - p.fraction)
        end
    end
        
    dxxdt_d[Nspecies*2] = receptorSurfpYCalc(xx_d, p.rr, AXL) + receptorSurfpYCalc(xx_d, p.rr, Mer) + receptorSurfpYCalc(xx_d, p.rr, Tyro)
    dxxdt_d[Nspecies*2] *= p.rr.pYinternalize*p.fraction      
end
