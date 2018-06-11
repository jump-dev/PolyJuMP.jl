# PolyJuMP Data
type Data
    # Delegates for polynomial constraints created
    delegates::Vector{ConstraintDelegate}
    # Default set for NonNegPoly
    nonnegpolydefault::Nullable
    # Default set for NonNegPolyMatrix
    nonnegpolymatrixdefault::Nullable
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
