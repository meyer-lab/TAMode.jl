using Test
using Profile
using TAMode
using LinearAlgebra
using ForwardDiff
using Turing

tps = 10.0 .^ range(-6.0, stop = 4.0, length = 25)
params = ones(15) * 0.5


aboutZero = x -> isapprox(x, 0.0, rtol = 1.0e-5, atol = 1.0e-5)


@testset "Can successfully assemble the parameters." begin
    TAMode.runTAM(tps, params, 10.0)

    @time TAMode.runTAM(tps, params, 10.0)
end


function gradFunc(x)
	return dot(TAMode.runTAM([10.0], x, 1.0), TAMode.pY)
end


@testset "Check that we can diff though the solver." begin
	g = similar(params)

    ForwardDiff.gradient!(g, gradFunc, params)
    @time ForwardDiff.gradient!(g, gradFunc, params)

    @test(length(g) == length(params))
end


include("testreactCode.jl")
include("testCompModel.jl")
include("testLS.jl")
include("testSample.jl")
