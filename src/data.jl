# PolyJuMP Data
mutable struct Data
    # Delegates for polynomial constraints created
    delegates::Vector{ConstraintDelegate}
    # Default set for NonNegPoly
    nonnegpolydefault
    # Default set for NonNegPolyMatrix
    nonnegpolymatrixdefault
    function Data()
        new(ConstraintDelegate[], nothing, nothing)
    end
end

function getpolydata(m::JuMP.Model)
    if !haskey(m.ext, :Poly)
        m.ext[:Poly] = Data()
    end
    m.ext[:Poly]
end
getdelegates(m::JuMP.Model) = getpolydata(m).delegates
