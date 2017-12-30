setpolymodule!(m::JuMP.Model, pm::Module) = setpolymodule!(getpolydata(m), pm)
setpolymodule!(data::PolyData, pm::Module) = data.polymodule = pm

getpolymodule(m::JuMP.Model) = getpolymodule(getpolydata(m))
function getpolymodule(data::PolyData)
    if isnull(data.polymodule)
        return DefaultModule
    end
    get(data.polymodule)
end

module DefaultModule
using JuMP, PolyJuMP, MultivariatePolynomials, SemialgebraicSets
export createpoly, polytype, addpolyconstraint!

polytype(m::JuMP.Model, p) = _polytype(m, p, p.x)
polytype(m::JuMP.Model, p, X::AbstractVector) = _polytype(m, p, monovec(X))

# Free polynomial

_polytype(m::JuMP.Model, ::Poly{false}, x::AbstractVector{MT}) where MT<:AbstractMonomial = polynomialtype(MT, JuMP.Variable)

# x should be sorted and without duplicates
function _createpoly(m::JuMP.Model, ::Poly{false}, x::AbstractVector{<:AbstractMonomial}, category::Symbol)
    polynomial((i) -> Variable(m, -Inf, Inf, category), x)
end

function createpoly(m::JuMP.Model, p::Poly{false, :Gram}, category::Symbol)
    _createpoly(m, p, monomials(sum(p.x)^2), category)
end

function createpoly(m::JuMP.Model, p::Union{Poly{false, :Default,Array{T}}, Poly{false, :Classic,Array{T}}}, category::Symbol) where T <: MultivariatePolynomials.AbstractMonomial
    _createpoly(m, p, monovec(p.x), category)
end

function createpoly(m::JuMP.Model, p::Union{PolyJuMP.Poly{false, :Default}, PolyJuMP.Poly{false, :Classic}}, category::Symbol)
    MultivariatePolynomials.polynomial([JuMP.Variable(m, -Inf, Inf, category)*p.x[i] for i in 1:length(p.x)])
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::FullSpace)
    constraints = JuMP.constructconstraint!.(coefficients(p), :(==))
    JuMP.addVectorizedConstraint(m, constraints)
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::AbstractAlgebraicSet)
    addpolyconstraint!(m, rem(p, ideal(domain)), s, FullSpace())
end

function addpolyconstraint!(m::JuMP.Model, p, s::ZeroPoly, domain::BasicSemialgebraicSet)
    addpolyconstraint!(m,  p, NonNegPoly(), domain)
    addpolyconstraint!(m, -p, NonNegPoly(), domain)
    nothing
end

addpolyconstraint!(m::JuMP.Model, args...) = error("PolyJuMP is just a JuMP extension for modelling Polynomial Optimization: it does not implement any reformulation. To use automatic sums of squares (SOS) reformulations, install the SumOfSquares Julia package and try \`using SumOfSquares\` and \`setpolymodule!(SumOfSquares)\` or use \`SOSModel\` instead of \`Model\`.")

end
