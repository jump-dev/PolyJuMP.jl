module TestFunctions

using Test
import MutableArithmetics as MA

using MultivariatePolynomials
const MP = MultivariatePolynomials

using JuMP
using PolyJuMP

function test_functions(var)
    x = MP.similar_variable(var, Val{:x})
    α = MOI.VariableIndex(1)
    @test α * x isa AbstractTerm{MOI.VariableIndex}
    @test x * α isa AbstractTerm{MOI.VariableIndex}
    @test α * x^2 isa AbstractTerm{MOI.VariableIndex}
    @test x^2 * α isa AbstractTerm{MOI.VariableIndex}
    @test (1α) * x^2 isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test (1x)^2 * α isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test α * (1x) isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test (1x) * α isa AbstractTerm{MOI.ScalarAffineFunction{Int}}
    @test α * (1.0x) isa AbstractTerm{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x) * α isa AbstractTerm{MOI.ScalarAffineFunction{Float64}}
    @test α * (1x + 2) isa AbstractPolynomial{MOI.ScalarAffineFunction{Int}}
    @test (1x + 2) * α isa AbstractPolynomial{MOI.ScalarAffineFunction{Int}}
    @test α * (1.0x + 2.0) isa
          AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x + 2.0) * α isa
          AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
end

function test_coefficient_products(var)
    x = MP.similar_variable(var, Val{:x})
    α, β = MOI.VariableIndex(1), MOI.VariableIndex(2)
    for a in (α, 2α + 1), b in (1β, 3β + 2), p in (x, x + 1)
        q = b * p
        for product in (a * q, q * a)
            @test MA.isequal_canonical(MP.coefficient(product, x), a * b)
        end
    end
    return
end

function test_coefficient_sums(var)
    x = MP.similar_variable(var, Val{:x})
    α = MOI.VariableIndex(1)
    # MOI affine functions need matching numeric coefficient types.
    for a in (α, 2α + 1),
        p in (1x, 2x, x + 1, (2α + 1) * x, (2α + 1) * x + 1),
        op in (+, -)

        for value in (0, 2)
            @test MA.isequal_canonical(
                op(a, p)(x => value),
                op(a, p(x => value)),
            )
            @test MA.isequal_canonical(
                op(p, a)(x => value),
                op(p(x => value), a),
            )
        end
    end
    a = 1α
    for (result, op) in (
        (MA.@rewrite((2x + 3) + a * (x + 1)), +),
        (MA.@rewrite((2x + 3) - (x + 1) * a), -),
    )
        @test MA.isequal_canonical(MP.coefficient(result, x), op(2, a))
        @test MA.isequal_canonical(MP.coefficient(result, one(x)), op(3, a))
    end
    return
end

function test_bare_monomial_sums(var)
    x = MP.similar_variable(var, Val{:x})
    α = MOI.VariableIndex(1)
    # Evaluate MOI coefficients first: Base arithmetic does not support mixing
    # ScalarAffineFunction{Int} and Bool during polynomial evaluation.
    function evaluate(p)
        q = MP.map_coefficients(p) do c
            return MOI.Utilities.eval_variables(_ -> 2, c)
        end
        return q(x => 3)
    end
    @test evaluate(α + x) == 5
    @test evaluate(x + α) == 5
    # Raw MOI arithmetic infers Bool coefficients for bare monomials, which
    # cannot represent the negative coefficient of subtraction.
    @test_broken evaluate(α - x) == -1
    @test_broken evaluate(x - α) == 1
    # Rewriting also mixes Bool coefficients from α * x with Int coefficients.
    @test_broken evaluate(MA.@rewrite((2x + 3) + α * (x + 1))) == 17
    @test_broken evaluate(MA.@rewrite((2x + 3) - (x + 1) * α)) == 1
    return
end

function test_scalar_polynomial_function(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    z = MP.similar_variable(var, Val{:z})
    α = MOI.VariableIndex(1)
    β = MOI.VariableIndex(2)
    γ = MOI.VariableIndex(3)
    poly = x^2 + 2x * y + z^3
    f = PolyJuMP.ScalarPolynomialFunction{Int,typeof(poly)}(poly, [β, γ, α])
    g = MOI.Utilities.substitute_variables(f) do var
        i = var.value
        j = mod1(i + 1, 3)
        return MOI.VariableIndex(j) + i
    end
    a, b, c = MP.variables(g.polynomial)
    X = c + 2
    Y = a + 3
    Z = b + 1
    @test g.polynomial == X^2 + 2X * Y + Z^3
    @test g.variables == [α, β, γ]
    @test iszero(MOI.constant(f))
    @test MOI.constant(g) == 17
end

function runtests(var)
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name) $(typeof(var))" begin
                getfield(@__MODULE__, name)(var)
            end
        end
    end
end

end

import DynamicPolynomials
DynamicPolynomials.@polyvar(x)
TestFunctions.runtests(x)
import TypedPolynomials
TypedPolynomials.@polyvar(y)
TestFunctions.runtests(y)
