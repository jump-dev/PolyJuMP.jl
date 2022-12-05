using JuMP
const MOIT = MOI.Test

MOI.Utilities.@model(
    NoFreeVariable,
    (),
    (MOI.EqualTo, MOI.LessThan, MOI.GreaterThan),
    (MOI.Nonnegatives, MOI.Nonpositives, MOI.Zeros),
    (),
    (),
    (MOI.ScalarAffineFunction,),
    (MOI.VectorOfVariables,),
    (MOI.VectorAffineFunction,)
)
# No free variables to make sure variable bridges are used to increase coverage
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Reals},
)
    return false
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Nonnegatives},
)
    return true
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Nonpositives},
)
    return true
end
function MOI.supports_add_constrained_variables(
    ::NoFreeVariable,
    ::Type{MOI.Zeros},
)
    return true
end

function bridged_mock(
    mock_optimize!::Function...;
    model = MOI.Utilities.Model{Float64}(),
)
    mock = MOI.Utilities.MockOptimizer(model)
    bridged = MOI.Bridges.full_bridge_optimizer(mock, Float64)
    MOI.Utilities.set_mock_optimize!(mock, mock_optimize!...)
    return bridged
end
