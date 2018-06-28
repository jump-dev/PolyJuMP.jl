export setpolymodule!

setdefault!(m::JuMP.Model, S::Type, T::Type) = setdefault!(getpolydata(m), S, T)
setdefault!(data::Data, ::Type{NonNegPoly}, S::Type) = data.nonnegpolydefault = S
setdefault!(data::Data, ::Type{NonNegPolyMatrix}, S::Type) = data.nonnegpolymatrixdefault = S

getdefault(m::JuMP.Model, s) = getdefault(getpolydata(m), s)

# includes Poly but also custom types defined by other modules
getdefault(data::Data, p) = p
getdefault(data::Data, s::Union{NonNegPoly, NonNegPolyMatrix}) = getdefault(data, typeof(s))()
getdefault(data::Data, ::Type{NonNegPoly}) = get_default(data.nonnegpolydefault)
getdefault(data::Data, ::Type{NonNegPolyMatrix}) = get_default(data.nonnegpolymatrixdefault)

setpolymodule!(m::JuMP.Model, pm::Module) = setpolymodule!(getpolydata(m), pm)
setpolymodule!(data::Data, pm::Module) = pm.setdefaults!(data)

function get_default(d::Nullable)
    if isnull(d)
        error("PolyJuMP is just a JuMP extension for modelling Polynomial Optimization: it does not implement any reformulation. To use automatic sums of squares (SOS) reformulations, install the SumOfSquares Julia package and try \`using SumOfSquares\` and \`setpolymodule!(SumOfSquares)\` or use \`SOSModel\` instead of \`Model\`.")
    end
    get(d)
end
