@testset "Model" begin
    optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[-1.0], [-1.0]])
    mock = bridged_mock(optimize!)
    Tests.real_zero_polynomial_test(mock, MOI.Test.Config())
end

@testset "Model" begin
    optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0, 0.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[-1.0], [-1.0]])
    mock = bridged_mock(optimize!)
    Tests.complex_zero_polynomial_test(mock, MOI.Test.Config())
end

# The VectorOfVariables-in-ZeroPolynomialSet is bridged by VectorFunctionizeBridge
# since the free variable is bridged. This tests that the MomentsAttribute is
# passed by the VectorFunctionizeBridge.
@testset "NoFreeVariable" begin
    optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        (MOI.ScalarAffineFunction{Float64}, MOI.LessThan{Float64}) => [-1.0],
        (MOI.VectorAffineFunction{Float64}, MOI.Zeros) => [[-1.0], [-1.0]])
    mock = bridged_mock(optimize!, model = NoFreeVariable{Float64}())
    Tests.real_zero_polynomial_test(mock, MOI.Test.Config())
end
