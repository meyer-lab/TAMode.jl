using MCMCChains

@testset "Test Chains Package." begin
    chn = read("../chain-file-10_1.jls", Chains)
end