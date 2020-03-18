tpsA549 = @SVector Float64[60, 240]
gasA549 = @SVector Float64[64, 16, 4, 1, 0.25, 0]
pYA549 = @SMatrix [10.8 8.3; 7.4 7.1; 7.1 7.7; 4.6 8.2; 6.1 7.2; 7.5 7.5]
totA549 = @SMatrix [3443.11 3219.69; 3143.41 3353.82; 3018.88 3611.82; 2608.88 3448.21; 2690.24 3168.14; 2672.00 2672.00]
surfA549 = @SMatrix [0.206 0.239; 0.274 0.316; 0.281 0.251; 0.220 0.302; 0.256 0.281; 0.257 0.337]


@model AXLfit(pYDataExp, surfDataExp, totDataExp, tps, g6conc, ::Type{TV}=Vector{Float64}) where {TV} = begin    
    paramsA ~ MvLogNormal(fill(-6.0, 2), 0.01)
    paramsB ~ Truncated(LogNormal(-1.0, 0.01), 0.0, 1.0)
    paramsC ~ MvLogNormal(fill(-6.0, 12), 0.01)
    scale ~ LogNormal(-1.0, 0.1)

    params = vcat(paramsA, paramsB, paramsC)
    pYresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))
    totalresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))
    surfresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))

    # Setup reductions
    pYAXL = TAMode.pY .* TAMode.recpSpecific[1]
    surfAXL = TAMode.surface .* TAMode.recpSpecific[1]
    totAXL = TAMode.total .* TAMode.recpSpecific[1]

    for ii = 1:length(g6conc)
        data = TAMode.runTAM(tps, params, g6conc[ii])

        pYresids[:, ii] = (data * pYAXL)*scale
        surfresids[:, ii] = (data * surfAXL)*scale
        totalresids[:, ii] = (data * totAXL)*scale
    end

    pYresids = vec(pYresids .- transpose(pYDataExp))
    totalresids = vec(totalresids .- transpose(totDataExp))
    surfresids = vec(surfresids .- transpose(surfDataExp))

    muu = zeros(length(pYresids))
    stdd = ones(length(pYresids))
    pYresids ~ MvNormal(muu, stdd*std(pYresids))
    surfresids ~ MvNormal(muu, stdd*std(surfresids))
    totalresids ~ MvNormal(muu, stdd*std(totalresids))
end
