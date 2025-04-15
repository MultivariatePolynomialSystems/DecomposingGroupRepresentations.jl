@testset "Lie algebras" begin
    @testset "LieAlgebra" begin
        so3 = algebra(LieGroup("SO", 3))
        @test size(so3) == 3
        @test dim(so3) == 3
        @test rank(so3) == 1
        @test nweights(so3) == 3
    end

    @testset "ScalingLieAlgebra" begin
        T = algebra(ScalingLieGroup([1 2 3 4; -1 -2 -3 -4]))
        @test size(T) == 4
        @test dim(T) == 2
        @test rank(T) == 2
    end

    @testset "Predefined weight structure of Lie algebra" begin
        so3 = algebra(LieGroup("SO", 3))
        ws = weight_structure(so3)
        cartan = cartan_subalgebra(so3)
        for wsp in ws
            for (Jel, λ) in zip(cartan, weight(wsp))
                Jm = matrix(Jel)
                @test Jm isa Matrix{ComplexF64}
                @test size(Jm) == (3, 3)
                for wv in wsp
                    v = vector(wv)
                    @test v isa Vector{ComplexF64}
                    @test length(v) == 3
                    @test Jm * vector(wv) ≈ λ * vector(wv)
                end
            end
        end
    end
end