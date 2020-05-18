rr = TAMode.Lsparam(fill(0.2, 9))
tps = @SVector Float64[60, 240]
GasProp = @SVector[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
GasConc = @SVector Float64[64, 16, 4, 1, 0.25, 0]
PROSConc = @SVector Float64[64, 16, 4, 1, 0.25, 0]


"Plot pY vs. ligand proportion."
function plotLigProp(rr, tps, GasProp)
    pY = Array{Float64}(undef, length(tps), length(GasProp))
    
    for i = 1:length(GasProp)
        gas = GasProp[i]
        pros = 1.0 - gas
        totalconc = 70
        pYdata = TAMode.runTAM(tps, rr, (gas*totalconc, pros*totalconc))
        pY[:, i] = pYdata * TAMode.pYLS
    end

    pY1 = [[pY[1, :], pY[2, :]]]
    plot(GasProp, pY1, label = ["1 hr" "4 hr"], title = "Phosphorylated Receptor vs. Ligand Proportion", lw = 3)
    xlabel!("Proportion of Gas6")
    ylabel!("Phosphorylated Receptor (pY)")
end


"Plot Gas6 dose response."
function plotGasDose(rr, tps, GasConc)
    pY = Array{Float64}(undef, length(tps), length(GasConc))
    
    for i = 1:length(GasConc)
        pYdata = TAMode.runTAM(tps, rr, (GasConc[i], 0.0))
        pY[:, i] = pYdata * TAMode.pYLS
    end

    pY2g = [[pY[1, :], pY[2, :]]]
    plot(GasConc, pY2g, label = ["1 hr" "4 hr"], title = "Gas6 Dose Response", lw = 3)
    xlabel!("Gas6 Concentration")
    ylabel!("Phosphorylated Receptor (pY)")
end


"Plot PROS1 dose response."
function plotPROSDose(rr, tps, PROSConc)
    pY = Array{Float64}(undef, length(tps), length(GasConc))
    
    for i = 1:length(PROSConc)
        pYdata = TAMode.runTAM(tps, rr, (0.0, PROSConc[i]))
        pY[:, i] = pYdata * TAMode.pYLS
    end

    pY2p = [[pY[1, :], pY[2, :]]]
    plot(PROSConc, pY2p, label = ["1 hr" "4 hr"], title = "PROS1 Dose Response", lw = 3)
    xlabel!("PROS1 Concentration")
    ylabel!("Phosphorylated Receptor (pY)")
end


"Plot dimer formation at 1 hr and 4 hr."
function plotDimers(rr, tps, GasProp)
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
        data = TAMode.runTAM(tps, rr, (gas*totalconc, pros*totalconc))
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
