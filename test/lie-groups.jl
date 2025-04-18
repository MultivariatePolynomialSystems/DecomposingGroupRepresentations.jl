@testset "Lie group Constructors" begin
    @testset "LieGroup" begin
        SO3 = LieGroup("SO", 3)
        @test SO3 isa LieGroup{ComplexF64}
        @test algebra(SO3) isa LieAlgebra{ComplexF64, Weight{Int}}
    end

    @testset "ScalingLieGroup" begin
        T = ScalingLieGroup([1 2 3 4; -1 -2 -3 -4])
        @test T isa ScalingLieGroup{ComplexF64}
        @test algebra(T) isa ScalingLieAlgebra{ComplexF64}
    end
end