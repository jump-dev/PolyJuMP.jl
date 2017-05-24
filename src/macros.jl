using JuMP
import JuMP: getvalue, validmodel, addtoexpr_reorder
using Base.Meta

export @set, @polyvariable

function getvalue{C}(p::Polynomial{C, JuMP.Variable})
    Polynomial(map(getvalue, p.a), p.x)
end
function getvalue{C}(p::MatPolynomial{C, JuMP.Variable})
    MatPolynomial(map(getvalue, p.Q), p.x)
end

polyvariable_error(args, str) = error("In @polyvariable($(join(args, ","))): ", str)

macro polyvariable(args...)
    length(args) <= 1 &&
    polyvariable_error(args, "Expected model as first argument, then variable information.")
    m = esc(args[1])
    x = args[2]
    extra = vcat(args[3:end]...)

    gottype = false
    haslb = false
    hasub = false
    # Identify the variable bounds. Five (legal) possibilities are "x >= lb",
    # "x <= ub", "lb <= x <= ub", "x == val", or just plain "x"
    explicit_comparison = false
    if isexpr(x,:comparison) # two-sided
        polyvariable_error(args, "Polynomial variable declaration does not support the form lb <= ... <= ub. Use ... >= 0 and separate constraints instead.")
    elseif isexpr(x,:call)
        explicit_comparison = true
        if x.args[1] == :>= || x.args[1] == :≥
            # x >= lb
            var = x.args[2]
            @assert length(x.args) == 3
            lb = JuMP.esc_nonconstant(x.args[3])
            if lb != 0
                polyvariable_error(args, "Polynomial variable declaration does not support the form ... >= lb with nonzero lb.")
            end
            nonnegative = true
        elseif x.args[1] == :<= || x.args[1] == :≤
            # x <= ub
            # NB: May also be lb <= x, which we do not support
            #     We handle this later in the macro
            polyvariable_error(args, "Polynomial variable declaration does not support the form ... <= ub.")
        elseif x.args[1] == :(==)
            polyvariable_error(args, "Polynomial variable declaration does not support the form ... == value.")
        else
            # Its a comparsion, but not using <= ... <=
            polyvariable_error(args, "Unexpected syntax $(string(x)).")
        end
    else
        # No bounds provided - free variable
        # If it isn't, e.g. something odd like f(x), we'll handle later
        var = x
        nonnegative = false
    end

    # separate out keyword arguments
    if VERSION < v"0.6.0-dev.1934"
        kwargs = filter(ex-> isexpr(ex,:kw), extra)
        extra  = filter(ex->!isexpr(ex,:kw), extra)
    else
        kwargs = filter(ex-> isexpr(ex,:(=)), extra)
        extra  = filter(ex->!isexpr(ex,:(=)), extra)
    end

    quotvarname = quot(getname(var))
    escvarname  = esc(getname(var))

    monotype = :None

    # process keyword arguments
    t = quot(:Cont)
    gram = false
    for ex in kwargs
        kwarg = ex.args[1]
        if kwarg == :grammonomials
            if monotype != :None
                polyvariable_error("Monomials given twice")
            end
            monotype = :Gram
            x = esc(ex.args[2])
        elseif kwarg == :monomials
            if monotype != :None
                polyvariable_error("Monomials given twice")
            end
            monotype = :Classic
            x = esc(ex.args[2])
        elseif kwarg == :category
            t = JuMP.esc_nonconstant(ex.args[2])
        else
            polyvariable_error(args, "Unrecognized keyword argument $kwarg")
        end
    end

    # Determine variable type (if present).
    # Types: default is continuous (reals)
    if isempty(extra)
        if monotype == :None
            polyvariable_error(args, "Missing monomial vector")
        end
    elseif length(extra) > 1
        polyvariable_error(args, "Too many extra argument: only expected monomial vector")
    else
        if monotype != :None
            polyvariable_error(args, "Monomials given twice")
        end
        monotype = :Default
        x = esc(extra[1])
    end
    Z = gensym()

    @assert monotype != :None

    create_fun = nonnegative ? :createnonnegativepoly : :createpoly

    # This is ugly I know
    monotypeid = monotype == :Default ? 1 : (monotype == :Classic ? 2 : 3)
    monotype = gensym()

    if isa(var,Symbol)
        # Easy case - a single variable
        return JuMP.assert_validmodel(m, quote
            $monotype = [:Default, :Classic, :Gram][$monotypeid]
            $escvarname = getpolymodule($m).$create_fun($m, $monotype, $x, $t)
        end)
    else
        polyvariable_error(args, "Invalid syntax for variable name: $(string(var))")
    end
end

function JuMP.constructconstraint!(p::Polynomial, sense::Symbol)
    PolyConstraint(sense == :(<=) ? -p : p, sense != :(==))
end

function appendconstraints!(domaineqs, domainineqs, expr::Expr, _error)
    if isexpr(expr, :call)
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
    elseif isexpr(expr, :&&)
        map(t -> appendconstraints!(domaineqs, domainineqs, t, _error), expr.args)
    else
        polyconstraint_error("Invalid domain constraint specification $(string(expr))")
    end
    nothing
end

function builddomain(domaineqs, domainineqs)
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
    domainaffs, code
end

macro set(expr)
    domaineqs = []
    domainineqs = []
    appendconstraints!(domaineqs, domainineqs, expr, msg -> error("In @set($expr: ", msg))
    domainvar, domaincode = builddomain(domaineqs, domainineqs)
    quote
        $domaincode
        $domainvar
    end
end
