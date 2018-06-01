using JuMP
using PolyJuMP
using Base.Test

using SemialgebraicSets

using MultivariatePolynomials
using DynamicPolynomials
#using TypedPolynomials

module TestPolyModule
using JuMP
using PolyJuMP
using MultivariatePolynomials
struct TestPoly{MS, MT<:MultivariatePolynomials.AbstractMonomial, MVT<:AbstractVector{MT}}
    x::MVT
end
TestPoly{MS}(x::AbstractVector{MT}) where {MS, MT} = TestPoly{MS, MT, typeof(x)}(x)
struct TestPolyVar{MS}
    x
    category::Symbol
end
JuMP.variabletype(m::JuMP.Model, p::TestPoly{MS}) where MS = TestPolyVar{MS}
PolyJuMP.createpoly(m::JuMP.Model, p::TestPoly{MS}, category::Symbol) where MS = TestPolyVar{MS}(p.x, category)

struct TestNonNegConstraint end
struct TestNonNegMatrixConstraint end
struct TestConstraint <: PolyJuMP.ConstraintDelegate
    p
    set
    domain
    kwargs
end
PolyJuMP.addpolyconstraint!(m::JuMP.Model, p, s::Union{TestNonNegConstraint, TestNonNegMatrixConstraint}, domain; kwargs...) = TestConstraint(p, s, domain, kwargs)

function setdefaults!(data::PolyJuMP.Data)
    PolyJuMP.setdefault!(data, PolyJuMP.Poly{true}, TestPoly)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPoly, TestNonNegConstraint)
    PolyJuMP.setdefault!(data, PolyJuMP.NonNegPolyMatrix, TestNonNegMatrixConstraint)
end

end

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

if isdir(Pkg.dir("SumOfSquares"))
    #include("sumofsquares.jl")
end
