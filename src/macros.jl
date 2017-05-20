using JuMP
import JuMP: getvalue, validmodel, addtoexpr_reorder
using Base.Meta

export Poly, @set

function getvalue{C}(p::Polynomial{C, JuMP.Variable})
    Polynomial(map(getvalue, p.a), p.x)
end
function getvalue{C}(p::MatPolynomial{C, JuMP.Variable})
    MatPolynomial(map(getvalue, p.Q), p.x)
end

# X is a vector of monomials to be used to construct a polynomial variable
# if MT is Gram, X represents the monomials of the form X^T Q X
# if MT is Classic, it represents the monomials of the form a^T X
# if MT is Default, it depends on whether the polynomials is constructed as nonnegative or not:
# For a nonnegative polynomial, it corresponds to Gram, otherwise it corresponds to Classic.
immutable Poly{P, MT, MV}
    x::MV
end
(::Type{Poly{P, MT}}){P, MT, MV}(x::MV) = Poly{P, MT, MV}(x)
(::Type{Poly{P}}){P}(x) = Poly{P, :Default}(x)
(::Type{Poly})(x) = Poly{false}(x)

function JuMP.variabletype(m::Model, p::Poly{true})
    getpolymodule(m).nonnegativepolytype(m, p)
end
function JuMP.variabletype(m::Model, p::Poly{false})
    getpolymodule(m).polytype(m, p)
end
function cvarchecks(_error::Function, lowerbound::Number, upperbound::Number, start::Number; extra_kwargs...)
    for (kwarg, _) in extra_kwargs
        _error("Unrecognized keyword argument $kwarg")
    end
    if !isnan(start)
        _error("Polynomial variable declaration does not support the form ... == value.")
    end
    if lowerbound != -Inf && upperbound != Inf
        _error("Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
    end
    if lowerbound != -Inf && lowerbound != 0
        _error("Polynomial variable declaration does not support the form ... >= lb with nonzero lb.")
    end
    if upperbound != Inf
        _error("Polynomial variable declaration does not support the form ... <= ub.")
    end
end
function JuMP.constructvariable!(m::Model, p::Poly{true}, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; extra_kwargs...)
    cvarchecks(_error, lowerbound, upperbound, start; extra_kwargs...)
    getpolymodule(m).createnonnegativepoly(m, p, category == :Default ? :Cont : category)
end
function JuMP.constructvariable!(m::Model, p::Poly{false}, _error::Function, lowerbound::Number, upperbound::Number, category::Symbol, basename::AbstractString, start::Number; extra_kwargs...)
    cvarchecks(_error, lowerbound, upperbound, start; extra_kwargs...)
    if lowerbound != -Inf
        _error("Free polynomial variable declaration does not support the form ... >= 0, use SOSPoly(x) instead of Poly(x) to create Sum of Squares polynomials. Note that SOSPoly(x) creates the polynomial x^T Q x with Q symmetric positive semidefinite while Poly(x) creates the polynomial a^T x so the meaning of the vector of monomial x changes from Poly to SOSPoly.")
    end
    getpolymodule(m).createpoly(m, p, category == :Default ? :Cont : category)
end

function JuMP.constructconstraint!(p::Polynomial, sense::Symbol)
    PolyConstraint(sense == :(<=) ? -p : p, sense != :(==))
end

function appendconstraints!(domains, domaineqs, domainineqs, expr, _error)
    if isexpr(expr, :call)
        try
            sense, vectorized = JuMP._canonicalize_sense(expr.args[1])
            @assert !vectorized
            if sense == :(>=)
                push!(domainineqs, :($(expr.args[2]) - $(expr.args[3])))
            elseif sense == :(<=)
                push!(domainineqs, :($(expr.args[3]) - $(expr.args[2])))
            elseif sense == :(==)
                push!(domaineqs, :($(expr.args[2]) - $(expr.args[3])))
            else
                polyconstraint_error("Unrecognized sense $(string(sense)) in domain specification")
            end
        catch
            push!(domains, esc(expr))
        end
    elseif isexpr(expr, :&&)
        map(t -> appendconstraints!(domains, domaineqs, domainineqs, t, _error), expr.args)
    else
        push!(domains, esc(expr))
    end
    nothing
end

function builddomain(domains, domaineqs, domainineqs)
    domainaffs = gensym()
    code = :( $domainaffs = isempty($domainineqs) ? (isempty($domaineqs) ? FullSpace() : AlgebraicSet()) : BasicSemialgebraicSet() )
    for dom in domaineqs
        affname = gensym()
        newaffdomain, parsecodedomain = JuMP.parseExprToplevel(dom, affname)
        code = quote
            $code
            $affname = zero(Polynomial{true, Int})
            $parsecodedomain
            addequality!($domainaffs, $newaffdomain)
        end
    end
    for dom in domainineqs
        affname = gensym()
        newaffdomain, parsecodedomain = JuMP.parseExprToplevel(dom, affname)
        code = quote
            $code
            $affname = zero(Polynomial{true, Int})
            $parsecodedomain
            addinequality!($domainaffs, $newaffdomain)
        end
    end
    for dom in domains
        code = quote
            $code
            $domainaffs = $domainaffs ∩ $dom
        end
    end
    domainaffs, code
end

macro set(expr)
    domains = []
    domaineqs = []
    domainineqs = []
    appendconstraints!(domains, domaineqs, domainineqs, expr, msg -> error("In @set($expr: ", msg))
    domainvar, domaincode = builddomain(domains, domaineqs, domainineqs)
    quote
        $domaincode
        $domainvar
    end
end
