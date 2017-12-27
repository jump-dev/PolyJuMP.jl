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
import JuMP, PolyJuMP, MultivariatePolynomials
function createpoly(m::JuMP.Model, p::Union{PolyJuMP.Poly{false, :Default}, PolyJuMP.Poly{false, :Classic}}, category::Symbol)
    MultivariatePolynomials.polynomial([JuMP.Variable(m, -Inf, Inf, category)*p.x[i] for i in 1:length(p.x)])
end

addpolyconstraint!(m::JuMP.Model, args...) = error("PolyJuMP is just a JuMP extension for modelling Polynomial Optimization but it does not implements any reformulation. You might want to run \`Pkg.add(\"SumOfSquares\")\` to install the Sum of Squares reformulation. If it is installed you can do \`using SumOfSquares\` and then \`setpolymodule!(SumOfSquares)\` to use it or use \`SOSModel\` instead of \`Model\`.")
end
