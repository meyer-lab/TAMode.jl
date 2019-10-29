using Test
using Profile
using TAMode

println("Testing.")

tps = [0.1, 1.0, 10.0, 100.0]


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, ones(15) * 0.5)
    @time TAMode.runTAM(tps, ones(15) * 0.5)

    @profile TAMode.runTAM(tps, ones(15) * 0.5)

    Profile.print()
end
