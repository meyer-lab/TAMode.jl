using Pkg
Pkg.instantiate()
using TAMode
using Plots
using StaticArrays

rr = TAMode.Lsparam(fill(0.2, 9));
tps = @SVector Float64[60, 240];

GasProp = @SVector[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
pYprop = Array{Float64}(undef, length(tps), length(GasProp));

GasConc = @SVector Float64[64, 16, 4, 1, 0.25, 0];
PROSConc = @SVector Float64[64, 16, 4, 1, 0.25, 0];
pY = Array{Float64}(undef, length(tps), length(GasConc));


#Phosphorylated Receptor vs. Ligand Proportion
for i = 1:length(GasProp)
    gas = GasProp[i]
    pros = 1.0 - gas
    totalconc = 70
    pYdata = TAMode.runTAM(tps, rr, (gas*totalconc, pros*totalconc))
    pYprop[:, i] = pYdata * TAMode.pYLS
end

pY1 = [[pYprop[1, :], pYprop[2, :]]]
plot(GasProp, pY1, label = ["1 hr" "4 hr"], title = "Phosphorylated Receptor vs. Ligand Proportion", lw = 3)
xlabel!("Proportion of Gas6")
ylabel!("Phosphorylated Receptor (pY)")


#Gas6 Dose Response
for i = 1:length(GasConc)
    pYdata = TAMode.runTAM(tps, rr, (GasConc[i], 0.0))
    pY[:, i] = pYdata * TAMode.pYLS
end

pY2g = [[pY[1, :], pY[2, :]]]
plot(GasConc, pY2g, label = ["1 hr" "4 hr"], title = "Gas6 Dose Response", lw = 3)
xlabel!("Gas6 Concentration")
ylabel!("Phosphorylated Receptor (pY)")


#PROS1 Dose Response
for i = 1:length(PROSConc)
    pYdata = TAMode.runTAM(tps, rr, (0.0, PROSConc[i]))
    pY[:, i] = pYdata * TAMode.pYLS
end

pY2p = [[pY[1, :], pY[2, :]]]
plot(PROSConc, pY2p, label = ["1 hr" "4 hr"], title = "PROS1 Dose Response", lw = 3)
xlabel!("PROS1 Concentration")
ylabel!("Phosphorylated Receptor (pY)")

