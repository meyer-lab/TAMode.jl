@testset "Test sampling." begin
    samp = sample(TAMode.AXLfit(TAMode.pYA549, TAMode.surfA549, TAMode.totA549, TAMode.tpsA549, TAMode.gasA549), HMC(0.001, 3), 2)
end
