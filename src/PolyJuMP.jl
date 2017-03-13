__precompile__()

module PolyJuMP

using MultivariatePolynomials
using JuMP
import JuMP: getdual, addconstraint
export getslack, setpolymodule!

# Polynomial Constraint

type PolyConstraint <: JuMP.AbstractConstraint
    p # typically either be a polynomial or a Matrix of polynomials
    nonnegative::Bool
    polymodule::Module
    domain::AbstractSemialgebraicSet
    delegate::Nullable
    function PolyConstraint(p, nonnegative::Bool, polymodule::Module, domain::AbstractSemialgebraicSet)
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

setpolymodule!(m::JuMP.Model, pm::Module) = setpolymodule!(getpolydata(m), pm)
setpolymodule!(data::PolyData, pm::Module) = data.polymodule = pm

getpolymodule(m::JuMP.Model) = getpolymodule(getpolydata(m))
function getpolymodule(data::PolyData)
    if isnull(data.polymodule)
        error("PolyJuMP is just a JuMP extension for modelling Polynomial Optimization but it does not implements any reformulation. You might want to run \`Pkg.add(\"SumOfSquares\")\` to install the Sum of Squares reformulation. If it is installed you can do \`using SumOfSquares\` and then \`setpolymodule!(SumOfSquares)\` to use it or use \`SOSModel\` instead of \`Model\`.")
    end
    get(data.polymodule)
end

include("macros.jl")
include("solve.jl")

end # module
