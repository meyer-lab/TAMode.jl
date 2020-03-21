tpsA549 = @SVector Float64[60, 240]
gasA549 = @SVector Float64[64, 16, 4, 1, 0.25, 0]
pYA549 = @SMatrix [10.8 8.3; 7.4 7.1; 7.1 7.7; 4.6 8.2; 6.1 7.2; 7.5 7.5]
totA549 = @SMatrix [3443.1 3219.7; 3143.4 3353.8; 3018.9 3611.8; 2608.9 3448.2; 2690.2 3168.1; 2672.0 2672.0]
surfA549 = @SMatrix [0.206 0.239; 0.274 0.316; 0.281 0.251; 0.220 0.302; 0.256 0.281; 0.257 0.337]


@model AXLfit(pYDataExp, surfDataExp, totDataExp, tps, g6conc, ::Type{TV} = Vector{Float64}) where {TV} = begin
    paramsA ~ MvLogNormal(fill(-3.0, 2), 0.01) # internalize, pYinternalize
    paramsB ~ Truncated(LogNormal(-1.0, 0.01), 0.0, 1.0) # sortF
    paramsC ~ MvLogNormal(fill(-3.0, 5), 0.01) # kRec, kDeg, xFwd, gasCur, AXLexpr
    Ig2rev ~ LogNormal(-3.0, 0.1)
    scale ~ LogNormal(-1.0, 0.1)
    scaleSurf ~ LogNormal(-1.0, 0.1)

    params = vcat(paramsA, paramsB, paramsC, zeros(2), Ig2rev, ones(4))
    pYresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))
    totalresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))
    surfresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))

    # Setup reductions
    pYAXL = TAMode.pY .* TAMode.recpSpecific[1]
    surfAXL = TAMode.surface .* TAMode.recpSpecific[1]
    totAXL = TAMode.total .* TAMode.recpSpecific[1]

    for ii = 1:length(g6conc)
        data = TAMode.runTAM(tps, params, g6conc[ii])

        pYresids[:, ii] = (data * pYAXL) * scale
        surfresids[:, ii] = (data * surfAXL) * scaleSurf
        totalresids[:, ii] = (data * totAXL)
    end

    pYresids = vec(pYresids .- transpose(pYDataExp))
    totalresids = vec(totalresids .- transpose(totDataExp))
    surfresids = vec(surfresids .- transpose(surfDataExp))

    muu = zeros(length(pYresids))
    stdd = ones(length(pYresids))
    pYresids ~ MvNormal(muu, stdd * std(pYresids))
    surfresids ~ MvNormal(muu, stdd * std(surfresids))
    totalresids ~ MvNormal(muu, stdd * std(totalresids))
end


A549model = AXLfit(TAMode.pYA549, TAMode.surfA549, TAMode.totA549, TAMode.tpsA549, TAMode.gasA549)
