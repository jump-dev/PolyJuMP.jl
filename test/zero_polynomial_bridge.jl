include("zero_polynomial_test.jl")
include("utilities.jl")

@testset "ZeroPolynomialBridge" begin
    config = MOI.Test.TestConfig()
    optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[-1.0]])
    mock = bridged_mock(optimize!)
    zero_polynomial_test(mock, config)
    F = MOI.VectorAffineFunction{Float64}
    S = PolyJuMP.ZeroPolynomialSet{FullSpace,MonomialBasis,Monomial{true},
                                   MonomialVector{true}}
    @test MOI.get(mock, MOI.ListOfConstraints()) == [
        (MOI.SingleVariable, MOI.LessThan{Float64}), (F, S)]
    @test MOI.get(mock, MOI.NumberOfConstraints{F, MOI.Zeros}()) == 0
    ci = first(MOI.get(mock, MOI.ListOfConstraintIndices{F, S}()))
    test_delete_bridge(mock, ci, 2, ((F, MOI.Zeros, 0),))
end
