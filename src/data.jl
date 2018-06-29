# PolyJuMP Data
type Data
    # Default set for NonNegPoly
    nonnegpolydefault::Nullable
    # Default set for NonNegPolyMatrix
    nonnegpolymatrixdefault::Nullable
    function Data()
        new(nothing, nothing)
    end
end

function getpolydata(m::JuMP.Model)
    if !haskey(m.ext, :Poly)
        m.ext[:Poly] = Data()
    end
    m.ext[:Poly]
end
