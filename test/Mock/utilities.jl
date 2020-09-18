using JuMP
const MOIT = MOI.Test

MOIU.@model(NoFreeVariable,
            (), (MOI.EqualTo, MOI.LessThan, MOI.GreaterThan), (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros), (),
            (), (MOI.ScalarAffineFunction,), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,))
# No free variables to make sure variable bridges are used to increase coverage
MOI.supports_add_constrained_variables(::NoFreeVariable, ::Type{MOI.Reals}) = false

function bridged_mock(mock_optimize!::Function...;
                      model = MOIU.Model{Float64}())
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end
