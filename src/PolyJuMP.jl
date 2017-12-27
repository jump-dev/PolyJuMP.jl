__precompile__()

module PolyJuMP

using MultivariatePolynomials
using SemialgebraicSets
using JuMP
export getslack, setpolymodule!

# Polynomial Constraint

export ZeroPoly, NonNegPoly

struct ZeroPoly end
struct NonNegPoly end

mutable struct PolyConstraint <: JuMP.AbstractConstraint
    p # typically either be a polynomial or a Matrix of polynomials
    set
    kwargs::Vector{Any}
    polymodule::Nullable{Module}
    domain::AbstractSemialgebraicSet
    delegate::Nullable
    function PolyConstraint(p, s)
        new(p, s, Any[], nothing, FullSpace(), nothing)
    end
end
function setpolymodule!(c::PolyConstraint, pm::Module)
    c.polymodule = pm
end
getpolymodule(c::PolyConstraint) = get(c.polymodule)

const PolyConstraintRef = ConstraintRef{Model, PolyConstraint}

function JuMP.addconstraint(m::Model, c::PolyConstraint; domain::AbstractSemialgebraicSet=FullSpace(), kwargs...)
    setpolymodule!(c, getpolymodule(m))
    c.domain = domain
    c.kwargs = kwargs
    polyconstr = getpolyconstr(m)
    push!(polyconstr, c)
    m.internalModelLoaded = false
    PolyConstraintRef(m, length(polyconstr))
end

function getdelegate(c::PolyConstraintRef, s::Symbol)
    delegate = getpolyconstr(c.m)[c.idx].delegate
    if isnull(delegate)
        error("$(string(s)) value not defined for constraint with index $c. Check that the model was properly solved.")
    end
    get(delegate)
end

function getslack(c::PolyConstraintRef)
    getslack(getdelegate(c, :Slack))
end
function JuMP.getdual(c::PolyConstraintRef)
    getdual(getdelegate(c, :Dual))
end

# PolyJuMP Data

type PolyData
    polyconstr::Vector{PolyConstraint}
    polymodule::Nullable{Module}
    function PolyData()
        new(PolyConstraint[], nothing)
    end
end

function getpolydata(m::JuMP.Model)
    if !haskey(m.ext, :Poly)
        m.solvehook = solvehook
        m.ext[:Poly] = PolyData()
    end
    m.ext[:Poly]
end

function getpolyconstr(m::JuMP.Model)
    getpolydata(m).polyconstr
end

include("macros.jl")
include("solve.jl")
include("module.jl")

end # module
