import Plots: plot, xlabel!, ylabel!, title!


"Plot pY vs. ligand proportion."
function plotLigProp(rr, tps)
    GasProp = LinRange(0.0, 1.0, 10)
    pY = Array{Float64}(undef, length(tps), length(GasProp))

    for i = 1:length(GasProp)
        gas = GasProp[i]
        pros = 1.0 - gas
        totalconc = 70
        pYdata = TAMode.runTAM(tps, rr, (gas * totalconc, pros * totalconc))
        pY[:, i] = pYdata * TAMode.pYLS
    end

    pY1 = [[pY[1, :], pY[2, :]]]
    plot(GasProp, pY1, label = ["1 hr" "4 hr"], title = "Phosphorylated Receptor vs. Ligand Proportion", lw = 3)
    xlabel!("Proportion of Gas6")
    ylabel!("Phosphorylated Receptor (pY)")
end


"Plot Gas6 dose response."
function plotGasDose(rr, tps, conc; Gas = true)
    pY = Array{Float64}(undef, length(tps), length(conc))

    for i = 1:length(conc)
        if Gas
            stim = (conc[i], 0.0)
        else
            stim = (0.0, conc[i])
        end

        pY[:, i] = TAMode.runTAM(tps, rr, stim) * TAMode.pYLS
    end

    pY2g = [[pY[1, :], pY[2, :]]]
    plot(conc, pY2g, label = ["1 hr" "4 hr"], title = "Gas6 Dose Response", lw = 3)
    xlabel!("Gas6 Concentration")
    ylabel!("Phosphorylated Receptor (pY)")
end


"Plot dimer formation at 1 hr and 4 hr."
function plotDimers(rr, tps)
    GasProp = LinRange(0.0, 1.0, 10)
    GGdimers = vcat(zeros(9), ones(2), zeros(3), zeros(9), ones(2), zeros(5))
    PPdimers = vcat(zeros(11), ones(2), zeros(12), ones(2), zeros(3))
    GPdimers = vcat(zeros(13), 1, zeros(13), 1, zeros(2))
    GG = Array{Float64}(undef, length(tps), length(GasProp))
    PP = Array{Float64}(undef, length(tps), length(GasProp))
    GP = Array{Float64}(undef, length(tps), length(GasProp))

    for i = 1:length(GasProp)
        gas = GasProp[i]
        pros = 1.0 - gas
        totalconc = 70
        data = TAMode.runTAM(tps, rr, (gas * totalconc, pros * totalconc))
        GG[:, i] = data * GGdimers
        PP[:, i] = data * PPdimers
        GP[:, i] = data * GPdimers
    end

    plot1hr = [[GG[1, :], PP[1, :]], GP[1, :]]
    plot(GasProp, plot1hr, label = ["GG" "PP" "GP"], lw = 3)
    title!("Dimer Formation at 1 hr")
    xlabel!("Proportion of Gas6")
    ylabel!("Dimer Concentration")

    plot4hr = [[GG[2, :], PP[2, :]], GP[2, :]]
    plot(GasProp, plot4hr, label = ["GG" "PP" "GP"], lw = 3)
    title!("Dimer Formation at 4 hr")
    xlabel!("Proportion of Gas6")
    ylabel!("Dimer Concentration")
end

"Plot pY receptor expression"
function plotpYExpression(chn, AXLexpr, MerTKexpr, Tyro3expr, cellName)

    tps = range(0.0, 1440, length = 200)

    pYAXL = TAMode.pY .* TAMode.recpSpecific[1]
    pYTyro3 = TAMode.pY .* TAMode.recpSpecific[3]

    data = calcData(chn, AXLexpr, MerTKexpr, Tyro3expr, tps)
    pYA = (data * pYAXL)  # normal pYA
    pYT = (data * pYTyro3)  # normal pYT

    # AXL knockdown
    A_data = calcData(chn, 0.0, MerTKexpr, Tyro3expr, tps)
    A_pYA = (A_data * pYAXL) # AXLkd pYA
    A_pYT = (A_data * pYTyro3) # AXLkd pYT

    # Tyro3 knockdown
    T_data = calcData(chn, AXLexpr, MerTKexpr, 0.0, tps)
    T_pYA = (T_data * pYAXL) # Tyrokd pYA
    T_pYT = (T_data * pYTyro3) # Tyrokd pYT

    pyAXL = vcat(pYA[1], pYA[2], A_pYA[1], A_pYA[2], T_pYA[1], T_pYA[2])
    pyTyro3 = vcat(pYT[1], pYT[2], A_pYT[1], A_pYT[2], T_pYT[1], T_pYT[2])

    count = vcat(pyAXL, pyTyro3)
    treatment = repeat(["Normal", "Normal", "AXL Knockdown", "AXL Knockdown", "Tyro3 Knockdown", "Tyro3 Knockdown"], outer = 2)
    time = repeat(["1 hr", "8 hr"], outer = 6)
    receptorType = repeat(["Phosphyrlated AXL", "Phosphyrlated Tyro3"], inner = 6)

    D = DataFrame(Treatment = treatment, Time = time, ReceptorType = receptorType, Count = count)

    p2 = Gadfly.plot(
        D,
        x = :Time,
        y = :Count,
        color = :ReceptorType,
        xgroup = :Treatment,
        Geom.subplot_grid(Geom.bar(position = :stack)),
        Guide.xlabel("Treatment"),
        Guide.title(cellName),
    )

end

"Plot expression over time"
function plotTimeSeries(chn, AXLexpr, MerTKexpr, Tyro3expr)

    tps = range(0.0, 1440, length = 200)
    data = calcData(chn, AXLexpr, MerTKexpr, Tyro3expr, tps)

    pYAXL = TAMode.pY .* TAMode.recpSpecific[1]
    surfAXL = TAMode.surface .* TAMode.recpSpecific[1]
    totAXL = TAMode.total .* TAMode.recpSpecific[1]

    pYMerTK = TAMode.pY .* TAMode.recpSpecific[2]
    surfMerTK = TAMode.surface .* TAMode.recpSpecific[2]
    totMerTK = TAMode.total .* TAMode.recpSpecific[2]

    pYTyro3 = TAMode.pY .* TAMode.recpSpecific[3]
    surfTyro3 = TAMode.surface .* TAMode.recpSpecific[3]
    totTyro3 = TAMode.total .* TAMode.recpSpecific[3]

    pYA = (data * pYAXL)
    surfA = (data * surfAXL)
    totalA = (data * totAXL)

    pYM = (data * pYMerTK)
    surfM = (data * surfMerTK)
    totalM = (data * totMerTK)

    pYT = (data * pYTyro3)
    surfT = (data * surfTyro3)
    totalT = (data * totTyro3)

    A = hcat(pYA, surfA, totalA)
    M = hcat(pYM, surfM, totalM)
    T = hcat(pYT, surfT, totalT)
    plot(
        tps,
        [A, M, T],
        lw = 3,
        xlabel = "Time (min)",
        label = ["AXL" "AXL" "AXL" "MerTK" "MerTK" "MerTK" "Tyro3" "Tyro3" "Tyro3"],
        title = ["pY" "surf" "total"],
        layout = (1, 3),
        size = (950, 400),
    )
end

function calcData(chn, AXLexpr, MerTKexpr, Tyro3expr, tps)
    index = 201
    Ig2rev = get(chn, :Ig2rev)[1]

    x = get(chn, [:internalize, :pYinternalize, :sortF, :kRec, :kDeg, :xFwd, :gasCur, :AXLexpr])
    samp_params = hcat(x.internalize, x.pYinternalize, x.sortF, x.kRec, x.kDeg, x.xFwd, x.gasCur)
    params = vcat(samp_params[index, :], AXLexpr, MerTKexpr, Tyro3expr, Ig2rev[index], [1.0, 1.0, 1.8, 100.0])

    data = TAMode.runTAM(tps, params, 10)
    return data
end


"Plot Comp Model."
function plotComp(pp, tps)
    r = 1:100
    pY = TAMode.compTAM(tps, pp)
    cplot = Array{Float64}(undef, length(tps), 100)

    for rr = 1:100
        for t = 1:length(tps)
            pYdata = view(pY, t, :, rr)
            cplot[:, rr] .= dot(pYdata, TAMode.pYc)
        end
    end

    plotpY = cplot[1, :]
    plot(r, plotpY, title = "Compartmental Model pY", lw = 3)
    
    if length(tps) > 1
        for tt = 2:length(tps)
            plotpY = cplot[tt, :]
            plot!(r, plotpY, lw = 3)
        end
    end
    xlabel!("Radius")
    ylabel!("pY")
    
end
