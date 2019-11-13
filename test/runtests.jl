using Test
using Profile
using TAMode

tps = [0.1, 1.0, 10.0, 100.0]
params = ones(15) * 0.5


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 1.0)
    @time TAMode.runTAM(tps, params, 1.0)

    @profile TAMode.runTAM(tps, params, 1.0)

    Profile.print(noisefloor=2.0)
end


@testset "Check for detailed balance at steady-state." begin
    uLong = TAMode.runTAM([1000.0], params, 100.0)

    dnorm = TAMode.TAM_reacti_dnorm(zeros(55), uLong, params, 0.0)

    # TODO: Lower tolerance
    @test dnorm < 2.0
end

@testset "Check if there's no receptor" begin
   
    
    for i=2:3
        rr::Rates = TAM.param(paramR.getTrafP()) #param is moved from reactcode.hpp
        
        rr.TAMs[i].expression = 0.0 #I don't know if we can use TAMs here... still working on this line
        
        tt::TAM = TAM(rr) #TAM is a constructor
        
        data::TAMout = tt.calcStim(tps,10)
        
        for t=0:size(tps)
            @test ≈ (data.total[i+3*t],0,1e-6)
        end
    end
end