struct ZeroPolynomialInAlgebraicSetBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    Z<:SA.AbstractBasis,
    DT<:SS.AbstractAlgebraicSet,
    B<:MB.SubBasis{MB.Monomial},
    Q<:MB.QuotientBasis,
} <: MOI.Bridges.Constraint.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{
        F,
        PolyJuMP.ZeroPolynomialSet{SS.FullSpace,Z,Q},
    }
    domain::DT
    basis::B
end

# The polynomial is zero on the algebraic set if it is in its ideal. The
# coefficients in the quotient basis are zero if and only if it is the case.
# If the Gröbner basis of the ideal is not known, it is not computed and the
# polynomial is instead constrained to be in the ideal truncated at its
# maximal degree, see `MultivariateMoments.quotient_basis`.
function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B,Q}},
    model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction,
    s::PolyJuMP.ZeroPolynomialSet{<:SS.AbstractAlgebraicSet,Z,B},
) where {T,F,Z,DT,B,Q}
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`,
    # we need to call `similar`. This is critical since `T` is `Float64` when
    # used with JuMP and the coefficient type is often `Int` with
    # `FixedVariablesSet`.
    quotient = MM.quotient_basis(s.basis, MP.similar(s.domain, T))
    coeffs = SA.coeffs(MOI.Utilities.scalarize(f), s.basis, quotient)
    zero_constraint = MOI.add_constraint(
        model,
        MOI.Utilities.vectorize(coeffs),
        PolyJuMP.ZeroPolynomialSet(SS.FullSpace(), s.zero_basis, quotient),
    )
    return ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B,Q}(
        zero_constraint,
        s.domain,
        s.basis,
    )
end

function MOI.supports_constraint(
    ::Type{ZeroPolynomialInAlgebraicSetBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{
        <:PolyJuMP.ZeroPolynomialSet{
            <:SS.AbstractAlgebraicSet,
            <:SA.AbstractBasis,
            <:MB.SubBasis{MB.Monomial},
        },
    },
) where {T}
    return true
end

# `FullSpace` is an `AbstractAlgebraicSet` but it is handled by
# `ZeroPolynomialBridge`
function MOI.supports_constraint(
    ::Type{ZeroPolynomialInAlgebraicSetBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{
        <:PolyJuMP.ZeroPolynomialSet{
            SS.FullSpace,
            <:SA.AbstractBasis,
            <:MB.SubBasis{MB.Monomial},
        },
    },
) where {T}
    return false
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ZeroPolynomialInAlgebraicSetBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B,Q}},
) where {T,F,Z,DT,B,Q}
    return [(F, PolyJuMP.ZeroPolynomialSet{SS.FullSpace,Z,Q})]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ZeroPolynomialInAlgebraicSetBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:PolyJuMP.ZeroPolynomialSet{DT,Z,B}},
) where {T,Z,DT<:SS.AbstractAlgebraicSet,B}
    G = MOI.Utilities.promote_operation(-, T, F, F)
    Q = MA.promote_operation(MM.quotient_basis, B, MP.similar_type(DT, T))
    return ZeroPolynomialInAlgebraicSetBridge{T,G,Z,DT,B,Q}
end

# Attributes, Bridge acting as an model
function MOI.get(
    ::ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B,Q},
    ::MOI.NumberOfConstraints{F,PolyJuMP.ZeroPolynomialSet{SS.FullSpace,Z,Q}},
) where {T,F,Z,DT,B,Q}
    return 1
end
function MOI.get(
    b::ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B,Q},
    ::MOI.ListOfConstraintIndices{
        F,
        PolyJuMP.ZeroPolynomialSet{SS.FullSpace,Z,Q},
    },
) where {T,F,Z,DT,B,Q}
    return [b.zero_constraint]
end

# Indices
function MOI.delete(model::MOI.ModelLike, c::ZeroPolynomialInAlgebraicSetBridge)
    return MOI.delete(model, c.zero_constraint)
end

# Attributes, Bridge acting as a constraint

function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintSet,
    bridge::ZeroPolynomialInAlgebraicSetBridge,
)
    set = MOI.get(model, attr, bridge.zero_constraint)
    return PolyJuMP.ZeroPolynomialSet(
        bridge.domain,
        set.zero_basis,
        bridge.basis,
    )
end
# TODO ConstraintPrimal

# The constraint is `A(p) = 0` where `A` maps the coefficients of `p` in
# `bridge.basis` to its coefficients in the quotient basis. The dual is
# therefore given by the adjoint map `A*`, see
# slide 9 of https://www.youtube.com/watch?v=C8dHxJCUHYw.
function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::ZeroPolynomialInAlgebraicSetBridge,
)
    dual = MOI.get(model, attr, bridge.zero_constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.zero_constraint)
    return SA.adjoint_coeffs(dual, bridge.basis, set.basis)
end

function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::ZeroPolynomialInAlgebraicSetBridge,
)
    return MOI.get(model, attr, bridge.zero_constraint)
end
