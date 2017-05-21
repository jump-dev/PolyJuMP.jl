export @polyconstraint

polyconstraint_error(m, x, args, str) = error("In @polyconstraint($m, $x$(isempty(args) ? "" : ", " * join(args, ","))): ", str)

macro polyconstraint(m, x, args...)
    warn("@polyconstraint is deprecated, use @constraint instead.")
    if isa(x, Symbol)
        polyconstraint_error(m, x, args, "Incomplete constraint specification $x. Are you missing a comparison (<= or >=)?")
    end

    (x.head == :block) &&
    polyconstraint_error(m, x, args, "Code block passed as constraint.")
    isexpr(x,:call) && length(x.args) == 3 || polyconstraint_error(m, x, args, "constraints must be in one of the following forms:\n" *
                                                                   "       expr1 <= expr2\n" * "       expr1 >= expr2")

    domaineqs = []
    domainineqs = []
    hasdomain = false
    for arg in args
        if ((VERSION < v"0.6.0-dev.1934") && !isexpr(arg, :kw)) || (VERSION >= v"0.6.0-dev.1934" && !isexpr(arg, :(=)))
            polyconstraint_error(m, x, args, "Unrecognized extra argument $(string(arg))")
        end
        if arg.args[1] == :domain
            @assert length(arg.args) == 2
            hasdomain && polyconstraint_error(m, x, args, "Multiple domain keyword arguments")
            hasdomain = true
            appendconstraints!(domaineqs, domainineqs, arg.args[2], msg -> polyconstraint_error(m, x, args, msg))
        else
            polyconstraint_error(m, x, args, "Unrecognized keyword argument $(string(arg))")
        end
    end

    # Build the constraint
    # Simple comparison - move everything to the LHS
    sense = x.args[1]
    if sense == :⪰
        sense = :(>=)
    elseif sense == :⪯
        sense = :(<=)
    end
    sense,_ = JuMP._canonicalize_sense(sense)
    lhs = :()
    if sense == :(>=) || sense == :(==)
        lhs = :($(x.args[2]) - $(x.args[3]))
    elseif sense == :(<=)
        lhs = :($(x.args[3]) - $(x.args[2]))
    else
        polyconstraint_error(m, x, args, "Invalid sense $sense")
    end
    newaff, parsecode = JuMP.parseExprToplevel(lhs, :q)
    domainaffs, domaincode = builddomain(domaineqs, domainineqs)
    code = quote
        q = zero(AffExpr)
        $parsecode
        $domaincode
    end
    nonnegative = !(sense == :(==))
    # I escape m here so that is is not escaped in the error messages of polyconstraint_error
    m = esc(m)
    JuMP.assert_validmodel(m, quote
        $code
        addconstraint($m, PolyConstraint($newaff, $nonnegative), domain = $domainaffs)
    end)
end
