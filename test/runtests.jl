using Test
using TAMode

println("Testing.")



@testset "Can successfully assemble the parameters." begin
    rr = TAMode.param(ones(15) * 0.5)
end