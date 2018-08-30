using DPSAGridTools

using Test
using SparseArrays
using LinearAlgebra

@time begin

        nodes = 4

        # use examples from ETH script
        # Example 2.1
        branch1 = line(1, 2, 0.0062, 0.036, 0.0105im)
        branch1p = line(1, 2, 0.0062, 0.036, 0.0105im) # two lines in parallel
        # Example 2.3
        branch2 = trafo_twoway(2, 3, 0, 0.23, 1.030)
        # Example 3.3, assume power base of 100MW
        branch3 = trafo_twoway(3, 4, 0, 0.0997, 1exp(1im*deg2rad(30)))
        # shunt at node 2
        shunt1 = Shunt(2, 0.1im)

        @testset "Grid Elements" begin

                # Example 2.1
                @test isequal(
                        2*round(branch1.circuit_matrix[1, 1] + branch1.circuit_matrix[1, 2] |> imag; sigdigits=5),
                        0.0105
                        )
                @test isequal(
                        imag(-branch1.circuit_matrix[1, 2]) / 0.0105 |> ceil,
                        -2569
                        )
                # Example 2.3
                @test isequal(
                        round(-branch2.circuit_matrix[1, 2] |> imag; sigdigits=3),
                        -4.48
                        )
                @test isequal(
                        round(branch2.circuit_matrix[2, 1] + branch2.circuit_matrix[2, 2] |> imag; sigdigits=3),
                        0.13
                        )
                # Example 3.3, assume power base of 100MW
                pu = 100 ./ [138^2 138*230; 230*138 230^2]
                V = [138 * 0.882 230*0.989exp(1im*deg2rad(-16.6))]
                Ikm = (branch3.circuit_matrix .* pu) * V'
                Pkm = V' .* conj(Ikm)
                @test isequal(
                        round(Pkm[1]; sigdigits=3),
                        203 - 70.8im
                )

        end

        @testset "Admittance" begin

                Y = build_nodal_admittance(nodes - 1, [branch1, branch1p, branch2, shunt1])
                # @show Array(Y)
                @test issymmetric(Y)

                Y = build_nodal_admittance(nodes, [branch1, branch1p, branch2, branch3, shunt1])
                # @show Array(Y)
                @test ~issymmetric(Y)

                @test typeof(Y) <: AbstractArray

                passive = [false, true, true, true]
                Yr = kron_reduction(Y, passive)

                @test typeof(Yr) <: AbstractArray

                @test isapprox(Yr[1], 2.4751958045854394e-5 - 0.07885604055123707im)

        end

end
