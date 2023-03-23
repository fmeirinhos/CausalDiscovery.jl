using CausalDiscovery
using Test

@testset verbose = true "CausalDiscovery.jl" begin
    @testset "skeleton" begin
        include("graph-inference.jl")
    end
end
