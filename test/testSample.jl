
@testset "Test sampling." begin
    samp = sample(TAMode.AXLfit(tpoints), HMC(0.01, 5), 10)
end
