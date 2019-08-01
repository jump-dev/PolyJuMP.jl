export setpolymodule!

setdefault!(m::JuMP.Model, S::Type, T::Type) = setdefault!(getpolydata(m), S, T)
setdefault!(data::Data, ::Type{NonNegPoly}, S::Type) = data.nonnegpoly_default = S
setdefault!(data::Data, ::Type{PosDefPolyMatrix}, S::Type) = data.posdefpolymatrix_default = S

getdefault(m::JuMP.Model, s) = getdefault(getpolydata(m), s)

# includes Poly but also custom types defined by other modules
getdefault(data::Data, p) = p
getdefault(data::Data, s::Union{NonNegPoly, PosDefPolyMatrix}) = getdefault(data, typeof(s))()
getdefault(data::Data, ::Type{NonNegPoly}) = get_default(data.nonnegpoly_default)
getdefault(data::Data, ::Type{PosDefPolyMatrix}) = get_default(data.posdefpolymatrix_default)

setpolymodule!(m::JuMP.Model, pm::Module) = setpolymodule!(getpolydata(m), pm)
setpolymodule!(data::Data, pm::Module) = pm.setdefaults!(data)

get_default(d) = d
function get_default(::Nothing)
    error("PolyJuMP is just a JuMP extension for modelling Polynomial Optimization: it does not implement any reformulation. To use automatic sums of squares (SOS) reformulations, install the SumOfSquares Julia package and try \`using SumOfSquares\` and \`setpolymodule!(model, SumOfSquares)\` or use \`SOSModel\` instead of \`Model\`.")
end
