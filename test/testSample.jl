@testset "Test sampling." begin
    samp = sample(TAMode.A549model, HMC(0.001, 10), 4)
end
