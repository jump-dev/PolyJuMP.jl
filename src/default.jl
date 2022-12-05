export setpolymodule!

setdefault!(m::JuMP.Model, S::Type, T::Type) = setdefault!(getpolydata(m), S, T)
function setdefault!(data::Data, ::Type{NonNegPoly}, S::Type)
    return data.nonnegpoly_default = S
end
function setdefault!(data::Data, ::Type{PosDefPolyMatrix}, S::Type)
    return data.posdefpolymatrix_default = S
end

getdefault(m::JuMP.Model, s) = getdefault(getpolydata(m), s)

# includes Poly but also custom types defined by other modules
getdefault(data::Data, p) = p
function getdefault(data::Data, s::Union{NonNegPoly,PosDefPolyMatrix})
    return getdefault(data, typeof(s))()
end
function getdefault(data::Data, ::Type{NonNegPoly})
    return get_default(data.nonnegpoly_default)
end
function getdefault(data::Data, ::Type{PosDefPolyMatrix})
    return get_default(data.posdefpolymatrix_default)
end

setpolymodule!(m::JuMP.Model, pm::Module) = setpolymodule!(getpolydata(m), pm)
setpolymodule!(data::Data, pm::Module) = pm.setdefaults!(data)

get_default(d) = d
function get_default(::Nothing)
    return error(
        "PolyJuMP is just a JuMP extension for modelling Polynomial Optimization: it does not implement any reformulation. To use automatic sums of squares (SOS) reformulations, install the SumOfSquares Julia package and try \`using SumOfSquares\` and \`setpolymodule!(model, SumOfSquares)\` or use \`SOSModel\` instead of \`Model\`.",
    )
end
