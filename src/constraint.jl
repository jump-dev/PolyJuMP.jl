export moments

abstract type PolynomialSet end

_in(::MIME) = Sys.iswindows() ? "in" : "∈"
_in(::MIME"text/latex") = "\\in"

function JuMP.in_set_string(mime::MIME, set::PolynomialSet)
    return string(_in(mime), ' ', set)
end

struct ZeroPoly <: PolynomialSet end
struct NonNegPoly <: PolynomialSet end
struct PosDefPolyMatrix <: PolynomialSet end

function JuMP.function_string(::MIME"text/plain", p::Union{SA.AlgebraElement,MP.AbstractPolynomialLike})
    return sprint(show, MIME"text/plain"(), p)
end
function JuMP.function_string(mime::MIME"text/latex", p::Union{SA.AlgebraElement,MP.AbstractPolynomialLike})
    return SA.strip_LaTeX(mime, sprint(show, MIME"text/latex"(), p))
end

### Shapes for polynomial/moments primal-dual pair ###

# Inspired from `JuMP.dual_shape` docstring example
struct PolynomialShape{B<:SA.ExplicitBasis} <: JuMP.AbstractShape
    basis::B
end
function JuMP.reshape_vector(x::Vector, shape::PolynomialShape)
    return MB.algebra_element(x, shape.basis)
end
struct MomentsShape{B<:SA.ExplicitBasis} <: JuMP.AbstractShape
    basis::B
end
function JuMP.reshape_vector(x::Vector, shape::MomentsShape)
    return MM.moment_vector(x, shape.basis)
end
JuMP.dual_shape(shape::PolynomialShape) = MomentsShape(shape.basis)
JuMP.dual_shape(shape::MomentsShape) = PolynomialShape(shape.basis)

JuMP.reshape_set(::ZeroPolynomialSet, ::PolynomialShape) = ZeroPoly()
function JuMP.moi_set(
    ::ZeroPoly,
    b::MB.SubBasis{MB.Monomial,M};
    domain::SS.AbstractSemialgebraicSet = SS.FullSpace(),
    basis = MB.FullBasis{MB.Monomial,M}(),
) where {M}
    return ZeroPolynomialSet(domain, basis, b)
end

"""
    moments(cref::JuMP.ConstraintRef)

Return the [`MomentsAttribute`](@ref) of `cref`.
"""
function MM.moments(cref::JuMP.ConstraintRef)
    return MOI.get(cref.model, MomentsAttribute(), cref)
end

"""
    bridges(F::Type{<:MOI.AbstractFunction}, S::Type{<:MOI.AbstractSet})

Return a list of bridges that may be needed to bridge `F`-in-`S` constraints but
not the bridges that may be needed by constraints added by the bridges.
"""
bridges(::Type{<:MOI.AbstractFunction}, ::Type{<:MOI.AbstractSet}) =
    Tuple{Type,Type}[]

"""
    bridges(S::Type{<:MOI.AbstractSet})

Return a list of bridges that may be needed to bridge variables constrained
in `S` on creation but not the bridges that may be needed by constraints added
by the bridges.
"""
bridges(S::Type{<:MOI.AbstractSet}) = bridges(MOI.VectorOfVariables, S)

function bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ZeroPolynomialSet{SS.FullSpace}},
)
    return Tuple{Type,Type}[(
        Bridges.Constraint.ZeroPolynomialBridge,
        coefficient_type_or_float(F),
    )]
end

function bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:ZeroPolynomialSet{<:SS.AbstractAlgebraicSet}},
)
    return Tuple{Type,Type}[(
        Bridges.Constraint.ZeroPolynomialInAlgebraicSetBridge,
        coefficient_type_or_float(F),
    )]
end

function bridges(F::Type{<:MOI.AbstractVectorFunction}, ::Type{<:PlusMinusSet})
    return Tuple{Type,Type}[(
        Bridges.Constraint.PlusMinusBridge,
        coefficient_type_or_float(F),
    )]
end

"""
    bridgeable(c::JuMP.AbstractConstraint, S::Type{<:MOI.AbstractSet})

Wrap the constraint `c` in `JuMP.BridgeableConstraint`s that may be needed to
bridge variables constrained in `S` on creation.

    bridgeable(c::JuMP.AbstractConstraint, F::Type{<:MOI.AbstractFunction},
               S::Type{<:MOI.AbstractSet})

Wrap the constraint `c` in `JuMP.BridgeableConstraint`s that may be needed to
bridge `F`-in-`S` constraints.
"""
function bridgeable end

function _concrete(
    bridge_type::Type{<:MOI.Bridges.Variable.AbstractBridge},
    S::Type{<:MOI.AbstractSet},
)
    return MOI.Bridges.Variable.concrete_bridge_type(bridge_type, S)
end
function _concrete(
    bridge_type::Type{<:MOI.Bridges.Constraint.AbstractBridge},
    S::Type{<:MOI.AbstractSet},
)
    return MOI.Bridges.Constraint.concrete_bridge_type(
        bridge_type,
        MOI.VectorOfVariables,
        S,
    )
end
function bridgeable(c::JuMP.AbstractConstraint, S::Type{<:MOI.AbstractSet})
    bridge_types = bridges(S)
    for bridge_type in bridge_types
        BT, T = bridge_type
        c = BridgeableConstraint(c, BT; coefficient_type = T)
        concrete_bridge_type = _concrete(BT{T}, S)
        for (ST,) in
            MOI.Bridges.added_constrained_variable_types(concrete_bridge_type)
            c = bridgeable(c, ST)
        end
        for (FT, ST) in MOI.Bridges.added_constraint_types(concrete_bridge_type)
            c = bridgeable(c, FT, ST)
        end
    end
    return c
end

# TODO If `JuMP.value_type` is not `Float64`, we should
#      redo `bridgeable`, `bridges` etc... in `JuMP.model_convert`
#      the most important is that `MOI.Bridges.concrete_bridge_type`
#      works so we should stick to the same coefficient type if there is
#      one. If there is none, `Float64` should be fine.
coefficient_type_or_float(::Type{<:MOI.AbstractFunction}) = Float64
coefficient_type_or_float(::Type{<:MOI.Utilities.TypedLike{T}}) where {T} = T

function bridgeable(
    c::JuMP.AbstractConstraint,
    F::Type{<:MOI.AbstractFunction},
    S::Type{<:MOI.AbstractSet},
)
    bridge_types = bridges(F, S)
    for bridge_type in bridge_types
        BT, T = bridge_type
        c = BridgeableConstraint(c, BT, coefficient_type = T)
        concrete_bridge_type =
            MOI.Bridges.Constraint.concrete_bridge_type(BT{T}, F, S)
        for (ST,) in
            MOI.Bridges.added_constrained_variable_types(concrete_bridge_type)
            c = bridgeable(c, ST)
        end
        for (FT, ST) in MOI.Bridges.added_constraint_types(concrete_bridge_type)
            c = bridgeable(c, FT, ST)
        end
    end
    return c
end

### @constraint macro ###

_vec(v::Vector) = v
_vec(v::AbstractVector) = collect(v)

non_constant_type(::Type{<:Real}) = AffExpr
non_constant_type(::Type{<:Complex}) = GenericAffExpr{ComplexF64,VariableRef}
non_constant_type(::Type{T}) where {T<:AbstractJuMPScalar} = T
function non_constant(a::AbstractVector{T}) where {T}
    # `coefficients` may return a `MP.LazyMap` and the `convert`
    # also take care of `collect`ing into a `Vector`
    return convert(Vector{non_constant_type(T)}, a)
end
non_constant_coefficients(p) = non_constant(MP.coefficients(p))

## ZeroPoly
function JuMP.build_constraint(
    error_fn::Function,
    p::MP.AbstractPolynomialLike,
    s::ZeroPoly;
    domain::SS.AbstractSemialgebraicSet = SS.FullSpace(),
    kws...,
)
    coefs = non_constant_coefficients(p)
    monos = MP.monomials(p)
    if domain isa SS.BasicSemialgebraicSet
        # p(x) = 0 for all x in a basic semialgebraic set. We replace it by
        # p(x) ≤ 0 and p(x) ≥ 0 for all x in the basic semialgebraic set.
        # We need to determine the cone two use for `NonNegPoly` which is stored
        # in the `ext` field of the `model` so we need the `model` hence we need
        # to go to `JuMP.add_constraint`.

        # We need to recreate the full list of keyword arguments. `kws` is a
        # `Iterators.Pairs, `values(kws)` is a `NamedTuple` and `keys(kws.itr)` is a
        # `Tuple`.
        all_kws = Iterators.Pairs(
            merge((domain = domain,), values(kws)),
            (:domain, values(kws)...),
        )
        return Constraint(error_fn, p, s, all_kws)
    else
        basis = MB.SubBasis{MB.Monomial}(monos)
        set = JuMP.moi_set(s, basis; domain = domain, kws...)
        constraint = JuMP.VectorConstraint(coefs, set, PolynomialShape(basis))
        return bridgeable(
            constraint,
            JuMP.moi_function_type(typeof(coefs)),
            typeof(set),
        )
    end
end
function JuMP.build_constraint(
    error_fn::Function,
    p::MP.AbstractPolynomialLike,
    s::MOI.EqualTo;
    kws...,
)
    return JuMP.build_constraint(error_fn, p - s.value, ZeroPoly(); kws...)
end

# # `NonNegPoly` and `PosDefPolyMatrix`
# The `model` is not given in `JuMP.build_constraint` so we create a custom
# `Constraint` object and transform the `set` in `JuMP.add_constraint`.
struct Constraint{F,PT,ST<:PolynomialSet,KWT<:Iterators.Pairs} <:
       JuMP.AbstractConstraint
    error_fn::F
    polynomial_or_matrix::PT
    set::ST
    kws::KWT
end
function JuMP.build_constraint(
    error_fn::Function,
    polynomial_or_matrix,
    set::Union{NonNegPoly,PosDefPolyMatrix};
    kws...,
)
    return Constraint(error_fn, polynomial_or_matrix, set, kws)
end

# FIXME the domain will not appear in the printing, it should be a field of
#       `ZeroPoly` which is `FullSpace` by default
JuMP.reshape_set(::PlusMinusSet, ::PolynomialShape) = ZeroPoly()
function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::Constraint{<:Any,<:Any,ZeroPoly},
    name::String = "",
)
    cone = getdefault(model, NonNegPoly())
    coefs = non_constant_coefficients(constraint.polynomial_or_matrix)
    monos = MP.monomials(constraint.polynomial_or_matrix)
    basis = MB.SubBasis{MB.Monomial}(monos)
    set = PlusMinusSet(JuMP.moi_set(cone, basis; constraint.kws...))
    new_constraint = JuMP.VectorConstraint(coefs, set, PolynomialShape(basis))
    bridgeable_con = bridgeable(
        new_constraint,
        JuMP.moi_function_type(typeof(coefs)),
        typeof(set),
    )
    return JuMP.add_constraint(model, bridgeable_con, name)
end

function JuMP.add_constraint(
    model::JuMP.Model,
    constraint::Constraint,
    name::String = "",
)
    set = getdefault(model, constraint.set)
    new_constraint = JuMP.build_constraint(
        constraint.error_fn,
        constraint.polynomial_or_matrix,
        set;
        constraint.kws...,
    )
    return JuMP.add_constraint(model, new_constraint, name)
end

# `NonNegPoly`
function JuMP.build_constraint(
    error_fn::Function,
    p::MP.AbstractPolynomialLike,
    s::MOI.GreaterThan;
    kws...,
)
    return JuMP.build_constraint(error_fn, p - s.lower, NonNegPoly(); kws...)
end
function JuMP.build_constraint(
    error_fn::Function,
    p::MP.AbstractPolynomialLike,
    s::MOI.LessThan;
    kws...,
)
    return JuMP.build_constraint(error_fn, s.upper - p, NonNegPoly(); kws...)
end

# `PosDefPolyMatrix`
# There is already a method for `AbstractMatrix` in `PSDCone` in `JuMP` so we
# need a more specific here to avoid ambiguity
function JuMP.build_constraint(
    error_fn::Function,
    p::AbstractMatrix{<:MP.AbstractPolynomialLike},
    s::PSDCone;
    kws...,
)
    return JuMP.build_constraint(error_fn, p, PosDefPolyMatrix(); kws...)
end

# Needed for the syntax `@constraint(model, A >= B, PSDCone())`
function JuMP.build_constraint(
    error_fn::Function,
    f::AbstractMatrix{<:MP.AbstractPolynomialLike},
    s::MOI.GreaterThan,
    extra::PSDCone,
)
    @assert iszero(s.lower)
    return JuMP.build_constraint(error_fn, f, extra)
end

# Needed for the syntax `@constraint(model, A <= B, PSDCone())`
function JuMP.build_constraint(
    error_fn::Function,
    f::AbstractMatrix{<:MP.AbstractPolynomialLike},
    s::MOI.LessThan,
    extra::PSDCone,
)
    @assert iszero(s.upper)
    new_f = MA.operate!!(*, -1, f)
    return JuMP.build_constraint(error_fn, new_f, extra)
end
