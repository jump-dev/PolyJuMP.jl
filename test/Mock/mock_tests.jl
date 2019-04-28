include("../Tests/Tests.jl")
include("utilities.jl")

using Test, JuMP

@testset "ZeroPolynomial" begin
    include("zero_polynomial.jl")
end
@testset "ZeroPolynomialInAlgebraicSet" begin
    include("zero_polynomial_in_algebraic_set.jl")
end
@testset "PlusMinus" begin
    include("plus_minus.jl")
end
