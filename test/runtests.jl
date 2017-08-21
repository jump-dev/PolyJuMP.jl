using DynamicPolynomials
using JuMP
using PolyJuMP
using Base.Test

using SemialgebraicSets

module TestPolyModule
using JuMP
using PolyJuMP
struct TestPoly{P}
    monotype::Symbol
    x
    category::Symbol
end
polytype{P, MT}(m::JuMP.Model, p::Poly{P, MT}) = TestPoly{P}
createpoly{P, MT}(m::JuMP.Model, p::Poly{P, MT}, category::Symbol) = TestPoly{P}(MT, p.x, category)
end

include("polymodule.jl")
include("variable.jl")
include("constraint.jl")

if isdir(Pkg.dir("SumOfSquares"))
    #include("sumofsquares.jl")
end
