module DummyPolyModule

using LinearAlgebra
import MathOptInterface as MOI
import StarAlgebras as SA
using MultivariatePolynomials
import MultivariateBases as MB
using SemialgebraicSets
using JuMP
using PolyJuMP

struct NonNeg{
    Z<:SA.AbstractBasis,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    B<:SA.ExplicitBasis,
} <: MOI.AbstractVectorSet
    zero_basis::Z
    domain::DT
    basis::B
    kwargs::Any
end
MOI.dimension(set::NonNeg) = length(set.monomials)
Base.copy(set::NonNeg) = set

struct DummyNonNeg <: PolyJuMP.PolynomialSet end

JuMP.reshape_set(::NonNeg, ::PolyJuMP.PolynomialShape) = DummyNonNeg()

struct DummyNonNegBridge{T,F} <: MOI.Bridges.Constraint.AbstractBridge end
function PolyJuMP.bridges(
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:NonNeg},
)
    return [(DummyNonNegBridge, PolyJuMP._coef_type(F))]
end
function MOI.Bridges.added_constrained_variable_types(
    ::Type{<:DummyNonNegBridge},
)
    return Tuple{Type}[]
end
function MOI.Bridges.added_constraint_types(
    ::Type{DummyNonNegBridge{T,F}},
) where {T,F}
    return Tuple{Type,Type}[(F, MOI.EqualTo{T})]
end
function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{<:DummyNonNegBridge{T}},
    F::Type{<:MOI.AbstractVectorFunction},
    ::Type{<:NonNeg},
) where {T}
    # This tests that `T` matches the coefficient type of `F`
    # as otherwise it would error.
    G = MOI.Utilities.promote_operation(-, T, F, MOI.VectorOfVariables)
    H = MOI.Utilities.scalar_type(G)
    return DummyNonNegBridge{T,H}
end

function JuMP.moi_set(
    cone::DummyNonNeg,
    b::MB.SubBasis{MB.Monomial,M};
    domain::AbstractSemialgebraicSet = FullSpace(),
    basis = MB.FullBasis{MB.Monomial,M}(),
    kwargs...,
) where {M}
    # For terms, `monomials` return a `OneOrZeroElementVector`
    # which would create a `NonNeg` of different type which would
    # make some complications in testing so let's convert it here
    return NonNeg(basis, domain, b, kwargs)
end

function JuMP.build_constraint(
    _error::Function,
    p::AbstractPolynomialLike,
    s::DummyNonNeg;
    kwargs...,
)
    coefs = PolyJuMP.non_constant_coefficients(p)
    basis = MB.SubBasis{MB.Monomial}(monomials(p))
    set = JuMP.moi_set(s, basis; kwargs...)
    return PolyJuMP.bridgeable(
        JuMP.VectorConstraint(coefs, set, PolyJuMP.PolynomialShape(basis)),
        JuMP.moi_function_type(typeof(coefs)),
        typeof(set),
    )
end

struct MatrixPolynomialShape{MT<:AbstractMonomial,MVT<:AbstractVector{MT}} <:
       JuMP.AbstractShape
    side_dimension::Int
    monomials::Matrix{MVT}
end

function JuMP.reshape_vector(
    x::Vector,
    shape::MatrixPolynomialShape{MT},
) where {MT}
    n = shape.side_dimension
    p = Matrix{polynomial_type(MT, eltype(x))}(undef, n, n)
    k = 0
    for j in 1:n
        for i in 1:n
            m = length(shape.monomials[i, j])
            p[i, j] = polynomial(x[k.+(1:m)], shape.monomials[i, j])
            k += m
        end
    end
    @assert length(x) == k
    return p
end

struct PosDefMatrix{
    Z<:SA.AbstractBasis,
    DT<:SemialgebraicSets.AbstractSemialgebraicSet,
    MT<:AbstractMonomial,
    MVT<:AbstractVector{MT},
} <: MOI.AbstractVectorSet
    zero_basis::Z
    domain::DT
    monomials::Matrix{MVT}
    kwargs::Any
end
MOI.dimension(set::PosDefMatrix) = sum(length, set.monomials)
Base.copy(set::PosDefMatrix) = set

struct DummyPosDefMatrix <: PolyJuMP.PolynomialSet end

function JuMP.reshape_set(::PosDefMatrix, ::MatrixPolynomialShape)
    return DummyPosDefMatrix()
end
function JuMP.moi_set(
    ::DummyPosDefMatrix,
    monos::Matrix{<:AbstractVector{<:AbstractMonomial}};
    domain::AbstractSemialgebraicSet = FullSpace(),
    basis = MB.FullBasis{MB.Monomial,eltype(eltype(monos))}(),
    kwargs...,
)
    return PosDefMatrix(basis, domain, monos, kwargs)
end

function JuMP.build_constraint(
    _error::Function,
    p::Matrix{<:AbstractPolynomialLike},
    s::DummyPosDefMatrix;
    kwargs...,
)
    n = LinearAlgebra.checksquare(p)
    # TODO we should use `non_constant_coefficients_type` once it exists
    coefs = coefficient_type(p[1, 1])[]
    for j in 1:n
        for i in 1:n
            append!(coefs, coefficients(p[i, j]))
        end
    end
    monos = monomials.(p)
    set = JuMP.moi_set(s, monos; kwargs...)
    return JuMP.VectorConstraint(coefs, set, MatrixPolynomialShape(n, monos))
end

# TODO the function in JuMP should not require the eltype to be
#      `AbstractJuMPScalar` so that we don't have to define this
# These methods are just copy-paste from JuMP/src/print.jl
function JuMP.function_string(
    ::MIME"text/plain",
    A::AbstractMatrix{<:AbstractPolynomialLike},
)
    str = sprint(show, MIME"text/plain"(), A)
    lines = split(str, '\n')
    # We drop the first line with the signature "m×n Array{...}:"
    lines = lines[2:end]
    # We replace the first space by an opening `[`
    lines[1] = '[' * lines[1][2:end]
    for i in 1:length(lines)
        lines[i] = lines[i] * (i == length(lines) ? ']' : ';')
    end
    return join(lines, '\n')
end
function JuMP.function_string(
    print_mode::MIME"text/latex",
    A::AbstractMatrix{<:AbstractPolynomialLike},
)
    str = sprint(show, MIME"text/plain"(), A)
    str = "\\begin{bmatrix}\n"
    for i in 1:size(A, 1)
        line = ""
        for j in 1:size(A, 2)
            if j != 1
                line *= " & "
            end
            if A isa Symmetric && i > j
                line *= "\\cdot"
            else
                line *= JuMP.function_string(print_mode, A[i, j])
            end
        end
        str *= line * "\\\\\n"
    end
    return str * "\\end{bmatrix}"
end

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, DummyNonNeg)
    return PolyJuMP.setdefault!(
        data,
        PolyJuMP.PosDefPolyMatrix,
        DummyPosDefMatrix,
    )
end

end
