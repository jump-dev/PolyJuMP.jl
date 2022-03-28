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
    @test α * (1.0x + 2.0) isa AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
    @test (1.0x + 2.0) * α isa AbstractPolynomial{MOI.ScalarAffineFunction{Float64}}
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
