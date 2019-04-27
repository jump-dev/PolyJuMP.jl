include("plus_minus_test.jl")
include("utilities.jl")
include("testpolymodule.jl")

@polyvar x

const NonNeg = TestPolyModule.NonNeg{MonomialBasis, typeof(@set x^2 ≤ 0),
                                     monomialtype(x), monovectype(x)}

MOIU.@model(PolyNonNegModel,
            (), (MOI.LessThan,), (NonNeg,), (),
            (MOI.SingleVariable,), (), (MOI.VectorOfVariables,),
            (MOI.VectorAffineFunction,))

@testset "PlusMinusBridge" begin
    config = MOI.Test.TestConfig()
    optimize!(mock) = MOIU.mock_optimize!(mock, [1.0, 1.0],
        (MOI.VectorOfVariables, NonNeg) => [[0.0]],
        (MOI.VectorAffineFunction{Float64}, NonNeg) => [[0.0, -0.5],
                                                        [0.0, 0.5],
                                                        [0.0]])
    mock = bridged_mock(optimize!; model = PolyNonNegModel{Float64}())
    plus_minus_test(mock, config; polymodule = TestPolyModule)
    F = MOI.VectorOfVariables
    G = MOI.VectorAffineFunction{Float64}
    S = PolyJuMP.PlusMinusSet{NonNeg}
    @test Set(MOI.get(mock, MOI.ListOfConstraints())) == Set([
        (MOI.SingleVariable, MOI.LessThan{Float64}), (G, S), (F, S)])
    ci = first(MOI.get(mock, MOI.ListOfConstraintIndices{F, S}()))
    test_delete_bridge(mock, ci, 2, ((F, NonNeg, 0), (G, NonNeg, 0));
                       last_bridge = false)
    ci = first(MOI.get(mock, MOI.ListOfConstraintIndices{G, S}()))
    test_delete_bridge(mock, ci, 2, ((F, NonNeg, 0), (G, NonNeg, 0)))
end
