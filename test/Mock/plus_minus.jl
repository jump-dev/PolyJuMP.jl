include("../testpolymodule.jl")

import MultivariateBases
const MB = MultivariateBases
using DynamicPolynomials
@polyvar x

using SemialgebraicSets
using PolyJuMP
const NonNeg = DummyPolyModule.NonNeg{
    MB.MonomialBasis,
    typeof(@set x^2 ≤ 0),
    monomial_type(x),
    monomial_vector_type(x),
}

MOI.Utilities.@model(
    PolyNonNegModel,
    (),
    (MOI.LessThan,),
    (NonNeg,),
    (),
    (),
    (),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,)
)

config = MOI.Test.Config()
function _optimize!(mock)
    return MOI.Utilities.mock_optimize!(
        mock,
        [1.0, 1.0],
        (MOI.VectorOfVariables, NonNeg) => [[0.0]],
        (MOI.VectorAffineFunction{Float64}, NonNeg) =>
            [[0.0, -0.5], [0.0, 0.5], [0.0]],
    )
end
mock = bridged_mock(_optimize!; model = PolyNonNegModel{Float64}())
Tests.plus_minus_test(mock, config; polymodule = DummyPolyModule)
