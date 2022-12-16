# PolyJuMP Data
mutable struct Data
    # Default set for NonNegPoly
    nonnegpoly_default::Any
    # Default set for PosDefPolyMatrix
    posdefpolymatrix_default::Any
    function Data()
        return new(nothing, nothing)
    end
end

function getpolydata(m::JuMP.Model)
    if !haskey(m.ext, :Poly)
        m.ext[:Poly] = Data()
    end
    return m.ext[:Poly]
end
