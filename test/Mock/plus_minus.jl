include("../testpolymodule.jl")

import MultivariatePolynomials as MP
import MultivariateBases as MB
using DynamicPolynomials
@polyvar x

using SemialgebraicSets
using PolyJuMP
const _full_basis_type = typeof(MB.FullBasis{MB.Monomial}(x))
const _sub_basis_type = MB.explicit_basis_type(_full_basis_type)
const NonNeg = DummyPolyModule.NonNeg{
    _full_basis_type,
    typeof(@set x^2 ≤ 0),
    _sub_basis_type,
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
            [[-0.5, 0.0], [0.5, 0.0], [0.0]],
    )
end
mock = bridged_mock(_optimize!; model = PolyNonNegModel{Float64}())
Tests.plus_minus_test(mock, config; polymodule = DummyPolyModule)
