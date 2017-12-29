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
