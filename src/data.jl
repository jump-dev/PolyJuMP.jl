# PolyJuMP Data
mutable struct Data
    # Default set for NonNegPoly
    nonnegpoly_default
    # Default set for PosDefPolyMatrix
    posdefpolymatrix_default
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
