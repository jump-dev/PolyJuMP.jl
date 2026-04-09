import MultivariatePolynomials as MP
import DynamicPolynomials
using IntervalArithmetic
const IA = IntervalArithmetic

const VariableOrder =
    DynamicPolynomials.Commutative{DynamicPolynomials.CreationOrder}
const MonomialOrder = MP.Graded{MP.LexOrder}
const VarType = DynamicPolynomials.Variable{VariableOrder,MonomialOrder}
const PolyType{T} = DynamicPolynomials.Polynomial{VariableOrder,MonomialOrder,T}
const FuncType{T} = PolyJuMP.ScalarPolynomialFunction{T,PolyType{T}}
const ExprMap = Dict{MOI.ScalarNonlinearFunction,MOI.VariableIndex}


function _to_reformulated_polynomial!(model::MOI.ModelLike, ::Type{T}, expr::Any) where {T}
	# println("creating ExprMap for expr: $expr")
	d = ExprMap()
	return _to_reformulated_polynomial!(model, T, d, expr)
end


function _to_reformulated_polynomial!(model::MOI.ModelLike, ::Type{T}, d::ExprMap, expr::Any) where {T}
	# println("Converting $expr, initially of type $(typeof(expr))")
	return _to_reformulated_polynomial!(model, T, d, convert(MOI.ScalarNonlinearFunction, expr))
end


function _to_reformulated_polynomial!(::MOI.ModelLike, ::Type{T}, ::ExprMap, expr::Union{MOI.VariableIndex,Number}) where {T}
	# println("Applying identity on $expr")
	return expr
end


function _to_reformulated_polynomial!(
    model::MOI.ModelLike,
    ::Type{T},
	d::ExprMap,
	expr::MOI.ScalarNonlinearFunction
) where {T}
	operator = expr.head
    operands = expr.args
	# println("I want to reformulate $expr")
    if operator == :*
		# in this version, you keep the first expr
		# new_operands = []
		# should_substitute = false
		# for operand in operands
		# 	if should_substitute
		# 		operand = _substitute_expr(model, T, d, operand)
		# 	end
		# 	push!(new_operands, operand)
		# 	if _is_expr(operand)
		# 		println("is an expr: $operand")
		# 		should_substitute = true
		# 	end
		# end
		# return MOI.ScalarNonlinearFunction(
		# 	:*,
		# 	new_operands
		# )

		# in this version, you substitude all expr products
		# if sum(_is_expr.(operands)) > 1
		# 	# println("$operands have more than one expr")
		# 	operands = _substitute_expr.(model, T, Ref(d), operands)
		# else
		# 	# println("$operands have less than two exprs")
		# 	operands = _to_reformulated_polynomial!.(model, T, Ref(d), operands)
		# end

		# in this version, operands are sorted by degree and then substituted from highest to lowest
		deg = _degree(expr)
		expr_count = sum(_is_expr.(operands))
		if deg > 2 && expr_count > 1
			sort!(operands, by=_degree, rev=true)
			for i in 1:length(operands)
				op = operands[i]
				# local_deg = _degree(op)
				# deg -= local_deg - 1
				operands[i] = _substitute_expr(model, T, d, op)
				expr_count -= 1
				if expr_count == 1
					break
				end
			end
		else
			operands = _to_reformulated_polynomial!.(model, T, Ref(d), operands)
		end
		return MOI.ScalarNonlinearFunction(
			:*,
			operands
		)
    elseif operator == :^ && length(operands) == 2
		a, b = operands
		# println("reformulation $a^$b...")
		# println("$a^$b is of degree $(_degree(expr))")
			#      if !(b isa Number && round(b) == b)
			#          b = round(b)
			# throw(error("EXPONENTIATION ERROR"))
			#      end
		deg = _degree(expr)
		num_vars = length(_var_index_set(expr))
		a = (deg < 3 || (deg == 3 && num_vars <= 2)) ? a : _substitute_expr(model, T, d, a)
		# a = (num_vars <= 1 || deg < 3 || (deg == 3 && num_vars <= 2)) ? a : _substitute_expr(model, T, d, a)
		# println("preparing to apply exponent: $a^$b, of types $(typeof(a)) and $(typeof(b))")
		return MOI.ScalarNonlinearFunction(:^, [a, b])
		# return a^b
    else
		return MOI.ScalarNonlinearFunction(
			operator,
			_to_reformulated_polynomial!.(model, T, Ref(d), operands)
		)
    end
end


_degree(expr::Number) = 0
_degree(expr::MOI.VariableIndex) = 1
_degree(expr::Any) = _degree(convert(MOI.ScalarNonlinearFunction, expr))

function _degree(expr::MOI.ScalarNonlinearFunction)
	operator = expr.head
	operands = expr.args

	if operator == :+ || operator == :-
		return maximum(_degree.(operands))
	elseif operator == :*
		return sum(_degree.(operands))
	elseif operator == :/ && length(operands) == 2
		return _degree(operands[1])
	elseif operator == :^ && length(operands) == 2
		if _degree(operands[2]) != 0
			throw(error("EXPONENTIATION ERROR"))
		end
		e = Int(_evaluate_func(operands[2]))
		return _degree(operands[1])*e
	else
		throw(error("Unknown expression $expr, cannot infer degree"))
	end
end


_evaluate_func(expr::MOI.VariableIndex) = throw(error("Cannot evaluate $expr, expr needs to be constant"))
_evaluate_func(expr::Number) = expr
_evaluate_func(expr::Any) = _evaluate_func(convert(MOI.ScalarNonlinearFunction, expr))

function _evaluate_func(expr::MOI.ScalarNonlinearFunction)
	operator = expr.head
	operands = expr.args
	if operator == :+
		return sum(_evaluate_func.(operands))
	elseif operator == :- && length(operands) == 2
		a, b = _evaluate_func.(operands)
		return a-b
	elseif operator == :- && length(operands) == 1
		return -_evaluate_func(operands[1])
	elseif operator == :*
		return prod(_evaluate_func.(operands))
	elseif operator == :/ && length(operands) == 2
		a, b = _evaluate_func.(operands)
		return a/b
	elseif operator == :^ && length(operands) == 2
		a, b = _evaluate_func.(operands)
		return a^b
	end
end


_var_index_set(::Number, ::Set{Int}) = nothing
_var_index_set(expr::MOI.VariableIndex, vars::Set{Int}) = push!(vars, expr.value)
_var_index_set(expr::Any, vars::Set{Int}) = _var_index_set(convert(MOI.ScalarNonlinearFunction, expr), vars)

function _var_index_set(expr::Any)::Set{Int}
	ret = Set{Int}()
	return _var_index_set(convert(MOI.ScalarNonlinearFunction, expr), ret)
end

function _var_index_set(expr::MOI.ScalarNonlinearFunction, vars::Set{Int})
	operands = expr.args
	for operand in operands
		_var_index_set(operand, vars)
	end
	return vars
end


function _substitute_expr(model::MOI.ModelLike, ::Type{T}, d::ExprMap, expr::Any) where {T}
	# println("sanity check")
	func = convert(MOI.ScalarNonlinearFunction, expr)
	if func.head ∉ [:+, :-]
		return _to_reformulated_polynomial!(model, T, d, func)
	end
	func = _to_reformulated_polynomial!(model, T, d, func)
	new_var = if haskey(d, func)
		# println("used cached variable $(d[func]) which is equal to $func")
		d[func]
	else
		bnds = _get_nonlinfunc_bounds(model, T, func)
		lb, ub = IA.bounds(bnds)
		# println("adding variable with bounds [$lb, $ub]")
		v, _ = MOI.add_constrained_variable(model, MOI.Interval(lb, ub))
		# v = MOI.add_variable(model)
		d[func] = v
		# println("added $func to dict")
		v
	end
	MOI.Utilities.normalize_and_add_constraint(model, MOI.ScalarNonlinearFunction(:-, [func, new_var]), MOI.EqualTo(0.0))
	# println("substitute $func by $new_var")
	return new_var
end


function _substitute_expr(::MOI.ModelLike, ::Type{T}, ::ExprMap, expr::Union{MOI.VariableIndex, Number}) where {T}
	return expr
end


function _is_expr(expr::MOI.ScalarNonlinearFunction)::Bool
	operator = expr.head
	operands = expr.args
	return (operator ∈ [:+, :-] && length(operands) > 1) || any(_is_expr, operands)
end

function _is_expr(::Union{MOI.VariableIndex,Number})::Bool
	return false
end

function _is_expr(expr::Any)::Bool
	return _is_expr(convert(MOI.ScalarNonlinearFunction, expr))
end


function _get_nonlinfunc_bounds(model::MOI.ModelLike, ::Type{T}, expr::Any)::IA.Interval{T} where T
	return _get_nonlinfunc_bounds(model, T, convert(MOI.ScalarNonlinearFunction, expr))
end


function _get_nonlinfunc_bounds(::MOI.ModelLike, ::Type{T}, expr::Number)::IA.Interval{T} where T
	return IA.interval(expr)
end


function _get_nonlinfunc_bounds(model::MOI.ModelLike, ::Type{T}, expr::MOI.VariableIndex)::IA.Interval{T} where T
	F = float(T)
	lb, ub = MOI.Utilities.get_bounds(model, T, expr)
	# println("got [$lb, $ub] for $expr")
	var_bounds = IntervalArithmetic.interval(
		lb == typemin(lb) ? typemin(F) : float(lb),
		ub == typemax(ub) ? typemax(F) : float(ub)
	)
	# println("returned $var_bounds")
	return var_bounds
end


function _get_nonlinfunc_bounds(model::MOI.ModelLike, ::Type{T}, expr::MOI.ScalarNonlinearFunction)::IA.Interval{T} where T
	operator = expr.head
	operands = expr.args
	if operator == :+
		return sum(_get_nonlinfunc_bounds.(model, T, operands))
	elseif operator == :-
		return -sum(_get_nonlinfunc_bounds.(model, T, operands))
	elseif operator == :*
		return prod(_get_nonlinfunc_bounds.(model, T, operands))
	elseif operator == :/
		return 1/prod(_get_nonlinfunc_bounds.(model, T, operands))
	elseif operator == :^ && length(operands) == 2
		a, b = _get_nonlinfunc_bounds.(model, T, operands)
		return a^b
	else
		throw(error("Unknown expression $expr"))
	end
end
