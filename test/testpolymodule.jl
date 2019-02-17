module TestPolyModule

using LinearAlgebra
using MathOptInterface
const MOI = MathOptInterface
using JuMP
using PolyJuMP
using MultivariatePolynomials
using SemialgebraicSets

struct NonNeg{BT <: PolyJuMP.AbstractPolynomialBasis,
              DT <: SemialgebraicSets.AbstractSemialgebraicSet,
              MT <: MultivariatePolynomials.AbstractMonomial,
              MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    domain::DT
    monomials::MVT
    kwargs
end
function Base.copy(set::NonNeg)
    return NonNeg(set.basis, set.domain, set.monomials, set.kwargs)
end

struct TestNonNeg <: PolyJuMP.PolynomialSet end

JuMP.reshape_set(::NonNeg, ::PolyJuMP.PolynomialShape) = TestNonNeg()
function JuMP.moi_set(cone::TestNonNeg,
                      monos::AbstractVector{<:AbstractMonomial};
                      domain::AbstractSemialgebraicSet=FullSpace(),
                      basis=MonomialBasis, kwargs...)
    return NonNeg(basis, domain, monos, kwargs)
end


function JuMP.build_constraint(_error::Function, p::AbstractPolynomialLike,
                               s::TestNonNeg; kwargs...)
    coefs = PolyJuMP.non_constant_coefficients(p)
    monos = monomials(p)
    set = JuMP.moi_set(s, monos; kwargs...)
    return JuMP.VectorConstraint(coefs, set, PolyJuMP.PolynomialShape(monos))
end

struct MatrixPolynomialShape{MT <: AbstractMonomial,
                             MVT <: AbstractVector{MT}} <: JuMP.AbstractShape
    side_dimension::Int
    monomials::Matrix{MVT}
end

function JuMP.reshape_vector(x::Vector,
                             shape::MatrixPolynomialShape{MT}) where {MT}
    n = shape.side_dimension
    p = Matrix{polynomialtype(MT, eltype(x))}(undef, n, n)
    k = 0
    for j in 1:n
        for i in 1:n
            m = length(shape.monomials[i, j])
            p[i, j] = polynomial(x[k .+ (1:m)], shape.monomials[i, j])
            k += m
        end
    end
    @assert length(x) == k
    return p
end

struct PosDefMatrix{BT <: PolyJuMP.AbstractPolynomialBasis,
                    DT <: SemialgebraicSets.AbstractSemialgebraicSet,
                    MT <: MultivariatePolynomials.AbstractMonomial,
                    MVT <: AbstractVector{MT}} <: MOI.AbstractVectorSet
    basis::Type{BT}
    domain::DT
    monomials::Matrix{MVT}
    kwargs
end
struct TestPosDefMatrix <: PolyJuMP.PolynomialSet end

function JuMP.reshape_set(::PosDefMatrix, ::MatrixPolynomialShape)
    return TestPosDefMatrix()
end
function JuMP.moi_set(::TestPosDefMatrix,
                      monos::Matrix{<:AbstractVector{<:AbstractMonomial}};
                      domain::AbstractSemialgebraicSet=FullSpace(),
                      basis=MonomialBasis, kwargs...)
    return PosDefMatrix(basis, domain, monos, kwargs)
end

function JuMP.build_constraint(_error::Function,
                               p::Matrix{<:AbstractPolynomialLike},
                               s::TestPosDefMatrix; kwargs...)
    n = LinearAlgebra.checksquare(p)
    # TODO we should use `non_constant_coefficients_type` once it exists
    coefs = coefficienttype(p[1, 1])[]
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
function JuMP.function_string(::Type{REPLMode},
                              A::AbstractMatrix{<:AbstractPolynomialLike})
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
function JuMP.function_string(print_mode::Type{IJuliaMode},
                              A::AbstractMatrix{<:AbstractPolynomialLike})
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
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, TestNonNeg)
    PolyJuMP.setdefault!(data, PolyJuMP.PosDefPolyMatrix, TestPosDefMatrix)
end

end
