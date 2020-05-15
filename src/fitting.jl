tpsA549 = @SVector Float64[60, 240]
gasA549 = @SVector Float64[64, 16, 4, 1, 0.25, 0]
pYA549 = @SMatrix [10.8 8.3; 7.4 7.1; 7.1 7.7; 4.6 8.2; 6.1 7.2; 7.5 7.5]
totA549 = @SMatrix [3443.1 3219.7; 3143.4 3353.8; 3018.9 3611.8; 2608.9 3448.2; 2690.2 3168.1; 2672.0 2672.0]
surfA549 = @SMatrix [0.206 0.239; 0.274 0.316; 0.281 0.251; 0.220 0.302; 0.256 0.281; 0.257 0.337]


" Run the model calculations for the A549 measurements. "
function dataModelCalc(tps, g6conc, params, scale, scaleSurf)
    pYresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))
    totalresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))
    surfresids = Matrix{typeof(scale)}(undef, length(tps), length(g6conc))

    # Setup reductions
    pYAXL = TAMode.pY .* TAMode.recpSpecific[1]
    surfAXL = TAMode.surface .* TAMode.recpSpecific[1]
    totAXL = TAMode.total .* TAMode.recpSpecific[1]

    params = param(params)
    solInit = getAutocrine(params)

    Threads.@threads for ii = 1:length(g6conc)
        params.gasCur = g6conc[ii]
        data = runTAMinit(tps, params, solInit)

        pYresids[:, ii] = (data * pYAXL) * scale
        surfresids[:, ii] = (data * surfAXL) * scaleSurf
        totalresids[:, ii] = (data * totAXL)
    end

    return pYresids, totalresids, surfresids
end


@model AXLfit(pYDataExp, surfDataExp, totDataExp, tps, g6conc, ::Type{TV} = Vector{Float64}) where {TV} = begin
    internalize ~ LogNormal(log(0.1), 0.1)
    pYinternalize ~ LogNormal(log(1.0), 0.1)
    sortF ~ Beta(1.0, 10.0)
    kRec ~ LogNormal(log(0.1), 0.1)
    kDeg ~ LogNormal(log(0.01), 0.1)
    xFwd ~ LogNormal(-3.0, 1.0)
    gasCur ~ LogNormal(log(0.1), 0.1)
    AXLexpr ~ LogNormal(log(1.0), 1.0)
    Ig2rev ~ LogNormal(0.0, 1.0)
    scale ~ LogNormal(0.0, 1.0)
    scaleSurf ~ LogNormal(0.0, 1.0)

    params = vcat(internalize, pYinternalize, sortF, kRec, kDeg, xFwd, gasCur, AXLexpr, zeros(2), Ig2rev, ones(4))
    pYresids, totalresids, surfresids = dataModelCalc(tps, g6conc, params, scale, scaleSurf)

    # The scaling is based on the average stderr, to make these proportional to the std MvNorm
    # TODO: Identify exact values
    sqResid = norm(pYresids .- transpose(pYDataExp))
    sqResid += norm(totalresids .- transpose(totDataExp)) / 100.0
    sqResid += norm(surfresids .- transpose(surfDataExp)) / 0.038

    sqResid ~ Chisq(length(pYDataExp) + length(totDataExp) + length(surfDataExp))
end


function plot_overlay(chn, tps, g6conc) 
    Ig2rev = get(chn, :Ig2rev)[1]
    scale = get(chn, :scale)[1]
    scaleSurf = get(chn, :scaleSurf)[1] 

    x = get(chn, [:internalize, :pYinternalize, :sortF, :kRec, :kDeg, :xFwd, :gasCur, :AXLexpr])
    samp_params = hcat(x.internalize, x.pYinternalize, x.sortF, x.kRec, x.kDeg, x.xFwd, x.gasCur, x.AXLexpr)

    pY = Array{Float64}(undef, size(samp_params, 1), length(tps), length(g6conc));
    tot = Array{Float64}(undef, size(samp_params, 1), length(tps), length(g6conc));
    surf = Array{Float64}(undef, size(samp_params, 1), length(tps), length(g6conc));
    
    for iter = 1:size(samp_params, 1)
        params = vcat(samp_params[iter, :], zeros(2), Ig2rev[iter], ones(4))

        pY[iter, :, :], tot[iter, :, :], surf[iter, :, :] = dataModelCalc(tps, g6conc, params, scale, scaleSurf)
    end

    # Calculate means
    medpY = Statistics.median(pY, dims = 1)
    medtot = Statistics.median(tot, dims = 1)
    medsurf = Statistics.median(surf, dims = 1)
    
    tp1_calcmed = hcat(transpose(medpY[:,1,:]), transpose(medsurf[:,1,:]), transpose(medtot[:,1,:]))
    tp2_calcmed = hcat(transpose(medpY[:,2,:]), transpose(medsurf[:,2,:]), transpose(medtot[:,2,:]))
    tp1_exp = hcat(pYA549[:,1], surfA549[:,1], totA549[:,1])
    tp2_exp = hcat(pYA549[:,2], surfA549[:,2], totA549[:,2])

    plot(gasA549, [tp1_calcmed, tp2_calcmed], 
        label=["1 hr, calc" "4 hr, calc"] , 
        title=["Phosphorylated receptor" "Surface receptor" "Total receptor"], 
        lw=3, 
        layout = (1,3), 
        size=(1200,400))
    plot!(gasA549, [tp1_exp, tp2_exp], 
        label=["1 hr, exp" "4 hr, exp"], 
        lw=3,
        layout=(1,3))
    xlabel!("Gas6 Concentration (nM)")
end

A549model = AXLfit(TAMode.pYA549, TAMode.surfA549, TAMode.totA549, TAMode.tpsA549, TAMode.gasA549)
