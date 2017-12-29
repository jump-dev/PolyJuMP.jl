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
