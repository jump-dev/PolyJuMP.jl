__precompile__()

module PolyJuMP

using MultivariatePolynomials
using JuMP
import JuMP: getdual, addconstraint
export getslack

# Polynomial Constraint

type PolyConstraint <: JuMP.AbstractConstraint
    p # typically either be a polynomial or a Matrix of polynomials
    nonnegative::Bool
    polymodule::Module
    domain::Vector
    delegate::Nullable
    function PolyConstraint(p, nonnegative::Bool, polymodule::Module, domain::Vector)
        new(p, nonnegative, polymodule, domain)
    end
end
PolyConstraint(p, nonnegative::Bool, polymodule::Module) = PolyConstraint(p, nonnegative, polymodule, Any[])

typealias PolyConstraintRef ConstraintRef{Model, PolyConstraint}

function addconstraint(m::Model, c::PolyConstraint)
    polyconstr = getpolyconstr(m)
    push!(polyconstr, c)
    m.internalModelLoaded = false
    PolyConstraintRef(m, length(polyconstr))
end

function getdelegate(c::PolyConstraintRef, s::Symbol)
    delegate = getpolyconstr(c.m)[c.idx].delegate
    if isnull(delegate)
        Base.warn("$(string(s)) value not defined for $(getname(v)). Check that the model was properly solved.")
    end
    get(delegate)
end

function getslack(c::PolyConstraintRef)
    getslack(getdelegate(c, :Slack))
end
function getdual(c::PolyConstraintRef)
    getdual(getdelegate(c, :Dual))
end

# PolyJuMP Data

type PolyData
    polyconstr::Vector{PolyConstraint}
    defaultpolymodule::Nullable{Module}
    function PolyData()
        new(PolyConstraint[], getdefaultpolymodule())
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

function getdefaultpolymodule(data::PolyData)
    if isnull(data.defaultpolymodule)
        error("PolyJuMP is just a JuMP extension for modelling Polynomial Optimization but it does not implements any reformulation. You might want to run \`Pkg.add(\"SumOfSquares\")\` to install the Sum of Squares reformulation.")
    end
    get(data.defaultpolymodule)
end

function getdefaultpolymodule(m::JuMP.Model)
    getdefaultpolymodule(getpolydata(m))
end

include("macros.jl")
include("solve.jl")
include("default.jl")

end # module
