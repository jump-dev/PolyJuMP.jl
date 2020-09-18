include("../testpolymodule.jl")

using DynamicPolynomials
@polyvar x

using SemialgebraicSets
using PolyJuMP
const NonNeg = TestPolyModule.NonNeg{MB.MonomialBasis, typeof(@set x^2 ≤ 0),
                                     monomialtype(x), monovectype(x)}

MOIU.@model(
    PolyNonNegModel,
    (), (MOI.LessThan,), (NonNeg,), (),
    (), (), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,))

config = MOI.Test.TestConfig()
_optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0],
    (MOI.VectorOfVariables, NonNeg) => [[0.0]],
    (MOI.VectorAffineFunction{Float64}, NonNeg) => [[0.0, -0.5],
                                                    [0.0, 0.5],
                                                    [0.0]])
mock = bridged_mock(_optimize!; model = PolyNonNegModel{Float64}())
Tests.plus_minus_test(mock, config; polymodule = TestPolyModule)
