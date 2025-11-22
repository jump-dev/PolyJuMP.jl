using LinearAlgebra

struct ZeroPolynomialInAlgebraicSetBridge{
    T,
    F<:MOI.AbstractVectorFunction,
    Z<:SA.AbstractBasis,
    DT<:SS.AbstractAlgebraicSet,
    B<:MB.SubBasis,
} <: MOI.Bridges.Constraint.AbstractBridge
    zero_constraint::MOI.ConstraintIndex{
        F,
        PolyJuMP.ZeroPolynomialSet{SS.FullSpace,Z,B},
    }
    domain::DT
    basis::B
end

function MOI.Bridges.Constraint.bridge_constraint(
    ::Type{ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B}},
    model::MOI.ModelLike,
    f::MOI.AbstractVectorFunction,
    s::PolyJuMP.ZeroPolynomialSet{
        <:SS.AbstractAlgebraicSet,
        Z,
        B,
    },
) where {T,F,Z,DT,B}
    p = MP.polynomial(MB.algebra_element(MOI.Utilities.scalarize(f), s.basis))
    # As `*(::MOI.ScalarAffineFunction{T}, ::S)` is only defined if `S == T`, we
    # need to call `similar`. This is critical since `T` is
    # `Float64` when used with JuMP and the coefficient type is often `Int` with
    # `FixedVariablesSet`.
    # FIXME convert needed because the coefficient type of `r` is `Any` otherwise if `domain` is `AlgebraicSet`
    r = convert(typeof(p), rem(p, SS.ideal(similar(s.domain, T))))
    zero_constraint = MOI.add_constraint(
        model,
        MOI.Utilities.vectorize(MP.coefficients(r)),
        PolyJuMP.ZeroPolynomialSet(
            SS.FullSpace(),
            s.zero_basis,
            MB.SubBasis{MB.Monomial}(MP.monomials(r)),
        ),
    )
    return ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B}(
        zero_constraint,
        s.domain,
        s.basis,
    )
end

function MOI.supports_constraint(
    ::Type{ZeroPolynomialInAlgebraicSetBridge{T}},
    ::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:PolyJuMP.ZeroPolynomialSet{<:SS.AbstractAlgebraicSet}},
) where {T}
    return true
end

function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:ZeroPolynomialInAlgebraicSetBridge},
)
    return Tuple{Type}[]
end

function MOI.Bridges.added_constraint_types(
    ::Type{<:ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B}},
) where {T,F,Z,DT,B}
    return [(
        F,
        PolyJuMP.ZeroPolynomialSet{
            SS.FullSpace,
            Z,
            B,
        },
    )]
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:ZeroPolynomialInAlgebraicSetBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:PolyJuMP.ZeroPolynomialSet{DT,Z,B}},
) where {T,Z,DT<:SS.AbstractAlgebraicSet,B}
    G = MOI.Utilities.promote_operation(-, T, F, F)
    return ZeroPolynomialInAlgebraicSetBridge{T,G,Z,DT,B}
end

# Attributes, Bridge acting as an model
function MOI.get(
    ::ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B},
    ::MOI.NumberOfConstraints{
        F,
        PolyJuMP.ZeroPolynomialSet{
            SS.FullSpace,
            Z,
            B,
        },
    },
) where {T,F,Z,DT,B}
    return 1
end
function MOI.get(
    b::ZeroPolynomialInAlgebraicSetBridge{T,F,Z,DT,B},
    ::MOI.ListOfConstraintIndices{
        F,
        PolyJuMP.ZeroPolynomialSet{
            SS.FullSpace,
            Z,
            B,
        },
    },
) where {T,F,Z,DT,B}
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

# Let A be the linear map corresponding to A(p) = rem(p, ideal(set.domain))
# We want to compute A*(μ) (where A* is the conjugate of A); see
# slide 9 of https://www.youtube.com/watch?v=C8dHxJCUHYw.
# That is, for every monomial of the basis `mono`, we want to compute
# `⟨mono, A*(μ)⟩`. By definition of the conjugacy, we have `⟨mono, A*(μ)⟩ = ⟨A(mono), μ⟩`
# so we can compute compute it by `⟨rem(mono, ideal(set.domain)), μ⟩`
function MOI.get(
    model::MOI.ModelLike,
    attr::MOI.ConstraintDual,
    bridge::ZeroPolynomialInAlgebraicSetBridge,
)
    dual = MOI.get(model, attr, bridge.zero_constraint)
    set = MOI.get(model, MOI.ConstraintSet(), bridge.zero_constraint)
    μ = MM.measure(dual, set.basis)
    I = SS.ideal(bridge.domain)
    return [dot(rem(mono, I), μ) for mono in MB.keys_as_monomials(bridge.basis)]
end

function MOI.get(
    model::MOI.ModelLike,
    attr::PolyJuMP.MomentsAttribute,
    bridge::ZeroPolynomialInAlgebraicSetBridge,
)
    return MOI.get(model, attr, bridge.zero_constraint)
end
