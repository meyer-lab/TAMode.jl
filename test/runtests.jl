using Test
using Profile
using TAMode
using LinearAlgebra
using ForwardDiff
using Turing

tps = 10.0 .^ range(-6.0, stop = 4.0, length = 25)

aboutZero = x -> isapprox(x, 0.0, rtol = 1.0e-5, atol = 1.0e-5)


include("testreactCode.jl")
include("testCompModel.jl")
include("testLS.jl")
include("testSample.jl")
