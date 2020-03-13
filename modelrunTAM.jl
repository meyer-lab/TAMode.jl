import TAMode
using Test
using Profile
using Turing

@model test(pYDataExp, surfDataExp, totDataExp, tps, g6conc, ::Type{TV}=Vector{Float64}) where {TV} = begin    
    paramsA ~ MvLogNormal(fill(-6.0, 2), 0.01)
    paramsB ~ Truncated(LogNormal(-1.0, 0.01), 0.0, 1.0)
    paramsC ~ MvLogNormal(fill(-6.0, 12), 0.01)
    scale ~ LogNormal(-1.0, 0.1)
    
    params = vcat(paramsA, paramsB, paramsC)
    resids = fill(Any[], 3, size(g6conc)[1])

    for index=1:size(g6conc)[1]
        #data = TAMode.getAutocrine(params)
        data = TAMode.runTAM(tps, params, g6conc[index]) #gas6 concentration

        #pY
        pYAXL = TAMode.pY .* TAMode.recpSpecific[1]
        pYData = data * pYAXL
        resids[1, index] = pYDataExp .- pYData*scale #scale?
        resids[1, index] ~ MvNormal(zeros(length(resids[1, index])), ones(length(resids[1, index]))*std(resids[1, index]))

        #surface
        surfAXL = TAMode.surface .* TAMode.recpSpecific[1]
        surfData = data * surfAXL
        resids[2, index] = surfDataExp .- surfData*scale
        resids[2, index] ~ MvNormal(zeros(length(resids[2, index])), ones(length(resids[2, index]))*std(resids[2, index]))

        #total
        totAXL = TAMode.total .* TAMode.recpSpecific[1]
        totData = data * totAXL
        resids[3, index] = totDataExp .- totData*scale
        resids[3, index] ~ MvNormal(zeros(length(resids[3, index])), ones(length(resids[3, index]))*std(resids[3, index]))
    end
end