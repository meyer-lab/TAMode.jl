tpsA549 = @SVector Float64[60, 240]
gasA549 = @SVector Float64[64, 16, 4, 1, 0.25, 0]
pYA549 = @SMatrix [10.8 8.3; 7.4 7.1; 7.1 7.7; 4.6 8.2; 6.1 7.2; 7.5 7.5]
totA549 = @SMatrix [3443.1 3219.7; 3143.4 3353.8; 3018.9 3611.8; 2608.9 3448.2; 2690.2 3168.1; 2672.0 2672.0]
surfA549 = @SMatrix [0.206 0.239; 0.274 0.316; 0.281 0.251; 0.220 0.302; 0.256 0.281; 0.257 0.337]


" Run the model calculations for the A549 measurements. "
function dataModelCalc(params, scale, scaleSurf)
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


" Full Turing A549 fitting model. "
@model AXLfit(pYDataExp, surfDataExp, totDataExp, tps, g6conc, ::Type{TV} = Vector{Float64}) where {TV} = begin
    internalize ~ LogNormal(-3.0, 0.01)
    pYinternalize ~ LogNormal(-3.0, 0.01)
    sortF ~ Truncated(LogNormal(-1.0, 0.01), 0.0, 1.0)
    kRec ~ LogNormal(-3.0, 0.01)
    kDeg ~ LogNormal(-3.0, 0.01)
    xFwd ~ LogNormal(-3.0, 0.01)
    gasCur ~ LogNormal(-3.0, 0.01)
    AXLexpr ~ LogNormal(-3.0, 0.01)
    Ig2rev ~ LogNormal(-3.0, 0.1)
    scale ~ LogNormal(-1.0, 0.1)
    scaleSurf ~ LogNormal(-1.0, 0.1)

    params = vcat(internalize, pYinternalize, sortF, kRec, kDeg, xFwd, gasCur, AXLexpr, zeros(2), Ig2rev, ones(4))
    pYresids, totalresids, surfresids = dataModelCalc(params, scale, scaleSurf)

    pYresids = vec(pYresids .- transpose(pYDataExp))
    totalresids = vec(totalresids .- transpose(totDataExp))
    surfresids = vec(surfresids .- transpose(surfDataExp))

    muu = zeros(length(pYresids))
    stdd = ones(length(pYresids))
    pYresids ~ MvNormal(muu, stdd * std(pYresids))
    surfresids ~ MvNormal(muu, stdd * std(surfresids))
    totalresids ~ MvNormal(muu, stdd * std(totalresids))
end


function plot_overlay(chn, tps, g6conc) 
    internalize = get(chn, :internalize)[1]
    pYinternalize = get(chn, :pYinternalize)[1]
    sortF = get(chn, :sortF)[1]
    kRec = get(chn, :kRec)[1]
    kDeg = get(chn, :kDeg)[1]
    xFwd = get(chn, :xFwd)[1]
    gasCur = get(chn, :gasCur)[1]
    AXLexpr = get(chn, :AXLexpr)[1]
    Ig2rev = get(chn, :Ig2rev)[1]
    scale = get(chn, :scale)[1]
    scaleSurf = get(chn, :scaleSurf)[1] 
    
    samp_params = hcat(internalize, pYinternalize, sortF, kRec, kDeg, xFwd, gasCur, AXLexpr)
    
    pY = Array{Float64}(undef, length(samp_params)[1], length(tps), length(g6conc));
    tot = Array{Float64}(undef, length(samp_params)[1], length(tps), length(g6conc));
    surf = Array{Float64}(undef, length(samp_params)[1], length(tps), length(g6conc));
    
    for iter = 1:length(samp_params)[1] 
        params = vcat(samp_params[iter,:], zeros(2), Ig2rev[iter], ones(4))

        pY[iter, :, :], tot[iter, :, :], surf[iter, :, :] = dataModelCalc(params, scale, scaleSurf)
    end
        
    #calculate means
    meanpY = fill(zeros(length(g6conc)), length(tps))
    meantot = fill(zeros(length(g6conc)), length(tps))
    meansurf = fill(zeros(length(g6conc)), length(tps))
    for time in 1:length(tps)
        for gas in 1:length(g6conc)
            meanpY[time][gas] = mean(pY[:, time, gas])
            meantot[time][gas] = mean(tot[:, time, gas])
            meansurf[time][gas] = mean(surf[:, time, gas])
        end
    end
    
    pYData = [pYA549[:,1],pYA549[:,2]] 
    surfData = [surfA549[:,1],surfA549[:,2]] 
    totData = [totA549[:,1],totA549[:,2]] 
        
    plot(g6conc, [meanpY; meansurf; meantot], 
            label=["1 hr, calc" "4 hr, calc"] , 
            title=["Phosphorylated receptor" "Surface receptor" "Total receptor"], 
            lw=3, 
            layout = (1,3), 
            size=(1200,400))
    plot!(g6conc, [pYData;surfData;totData], 
            label=["1 hr, exp" "4 hr, exp"], 
            lw=3)
    xlabel!("Gas6 Concentration (nM)")
    
end

A549model = AXLfit(TAMode.pYA549, TAMode.surfA549, TAMode.totA549, TAMode.tpsA549, TAMode.gasA549)