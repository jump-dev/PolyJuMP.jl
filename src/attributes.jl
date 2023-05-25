"""
    MomentsAttribute(N)
    MomentsAttribute()

A constraint attribute for the vector of moments corresponding to the
constraint that a polynomial is zero in the full space in result `N`. If the
polynomial is constrained to be zero in an algebraic set, it is the moments
for the constraint once it is rewritten into an constraint on the full space.
If `N` is omitted, it is 1 by default.

## Examples

Consider the following program:
```julia
@variable(model, α)
@variable(model, β ≤ 1)

using DynamicPolynomials
@polyvar x y
cref = @constraint(model, α * x - β * y == 0, domain = @set x == y)
```
The constraint is equivalent to
```julia
@constraint(model, (α - β) * y == 0)
```
for which the dual is the 1-element vector with the moment of `y` of value `-1`.
This is the result of `moments(cref)`. However, the dual of `cref` obtained by
`dual(cref)` is the 2-elements vector with both the moments of `x` and `y`
of value `-1`.
"""
struct MomentsAttribute <: MOI.AbstractConstraintAttribute
    result_index::Int
end
MomentsAttribute() = MomentsAttribute(1)

function MOI.is_set_by_optimize(::MomentsAttribute)
    return true
end

# If a variable is bridged, the `VectorOfVariables`-in-`ZeroPolynomialSet` is
# bridged by `MOI.Bridges.Constraint.VectorFunctionizeBridge` and it has
# to pass the constraint to the PolyJuMP bridge.
function MOI.Bridges.Constraint.invariant_under_function_conversion(
    ::MomentsAttribute,
)
    return true
end

# This is type piracy but we tolerate it.
const ObjectWithoutIndex = Union{
    AbstractMonomial,
    AbstractTerm{<:MOI.Utilities.ObjectWithoutIndex},
    AbstractPolynomial{<:MOI.Utilities.ObjectWithoutIndex},
    MB.AbstractPolynomialBasis,
    MultivariateMoments.AbstractMomentMatrix{
        <:MOI.Utilities.ObjectWithoutIndex,
    },
    MultivariateMoments.AbstractMeasure{<:MOI.Utilities.ObjectWithoutIndex},
}
const ObjectOrTupleWithoutIndex =
    Union{ObjectWithoutIndex,Tuple{Vararg{ObjectWithoutIndex}}}
const ObjectOrTupleOrArrayWithoutIndex =
    Union{ObjectOrTupleWithoutIndex,AbstractArray{<:ObjectOrTupleWithoutIndex}}
MOI.Utilities.map_indices(::Function, x::ObjectOrTupleOrArrayWithoutIndex) = x
function MOI.Utilities.substitute_variables(
    ::Function,
    x::ObjectOrTupleOrArrayWithoutIndex,
)
    return x
end
