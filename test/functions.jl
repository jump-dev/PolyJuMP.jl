module TestFunctions

using Test

using MultivariatePolynomials
const MP = MultivariatePolynomials

using JuMP
using PolyJuMP

function test_functions(var)
    x = MP.similarvariable(var, Val{:x})
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

function test_scalar_polynomial_function(var)
    x = MP.similarvariable(var, Val{:x})
    y = MP.similarvariable(var, Val{:y})
    z = MP.similarvariable(var, Val{:z})
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
