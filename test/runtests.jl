using Test
using Profile
using TAMode
using LinearAlgebra
using ForwardDiff
using Turing

tps = 10.0 .^ range(-6.0, stop = 4.0, length = 25)


#include("testreactCode.jl")
include("testCompModel.jl")
#include("testLS.jl")
#include("testSample.jl")
