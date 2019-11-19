using Test
using Profile
using TAMode

tps = [0.1, 1.0, 10.0, 100.0, 1000.0]
params = ones(15) * 0.5


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 1.0)
    @time TAMode.runTAM(tps, params, 1.0)

    @profile TAMode.runTAM(tps, params, 1.0)
    @profile TAMode.runTAM(tps, params, 10.0)
    @profile TAMode.runTAM(tps, params, 100.0)
    @profile TAMode.runTAM(tps, params, 1000.0)

    Profile.print(noisefloor=5.0)
end


@testset "Check for detailed balance at steady-state." begin
    uLong = TAMode.runTAM([1000000.0], params, 100.0)

    dnorm = TAMode.TAM_reacti_dnorm(zeros(55), uLong, params, 0.0)

    # TODO: Lower tolerance
    # Probably have to turn off trafficking to get this to work.
    @test dnorm < 1.7
end


@testset "Test that if none of the ligand is expressed, we don't end up seeing any." begin
    rr = TAMode.param(params)
    rr.gasCur = 0.0

    data = TAMode.runTAM(tps, rr, 0.0)

    @test all(data * TAMode.pY .≈ 0)
    @test all(data * TAMode.boundLig .≈ 0)
end


@testset "Check if there's no receptor that we don't see any." begin
    for i in 1:3
        rr = TAMode.param(params)
        rr.TAMs[i].expression = 0.0

        data = TAMode.runTAM(tps, rr, 1.0)

        @test all(data * TAMode.recpSpecific[i] .≈ 0.0)
    end
end

@testset "Ensure that system reaches equilibrium." begin
    tt = TAMode.param(params)
    tt.gasCur = tt.autocrine
    dnorm =  #??? think this is calling dnorm of parameters which calls l2 norm
    #tests that the probability density = 0 which in the C++ is the l2 norm version of the paired doubles
    #l2 norm calculates the squared root of the sum of the squared vector values
    @test dnorm < 1.7
end

#=
    void testTAMequilibrium() {
        TAM tt = TAM(TAM::Param(paramR.getTrafP()));
        tt.p.gasCur = tt.p.autocrine;

        CPPUNIT_ASSERT_NO_THROW(tt.initSystem());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(tt.getDNorm().first, 0, 1.0E-3);
    }
=#
