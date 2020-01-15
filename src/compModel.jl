compT = Array{Float64, 1}(undef, 13)

mutable struct comprates{T} #do we want this to be T and what does {T} mean
    rr::Rates #do i need values for all the rates and does this want to be in the reactCode file or a new one
    fraction::T #Fraction of cell surface covered with PtdSer
    partIn::T #Partitioning rate into PtdSer regions
    gasPart::T #Partitioning of ligand into PtdSer region
end


function compParamm(compIn::compT)
    paramsIn::trafT
    size = size(paramsIn)[1]
    for i in 1:size
        compIn[i+3] = paramsIn[i]
    end
    
    ppparams::comprates
    ppparams.rr = TAM.Param(paramsIn)
    ppparams.fraction = compIn[1]
    ppparams.partIn = compIn[2]
    ppparams.gasPart = compIn[3]
    
    return ppparams
end