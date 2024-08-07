module TestVariable

using Test

import MutableArithmetics as MA
import StarAlgebras as SA

using MultivariatePolynomials
const MP = MultivariatePolynomials
import MultivariateBases as MB
using DynamicPolynomials
#using TypedPolynomials

using JuMP
using PolyJuMP

function test_variable_macro_with_Poly(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    X = [1, y, x^2, x, y^2]

    m = Model()
    @variable(m, a)
    var_poly = polynomial(a * x * y)
    var_poly_x = polynomial(a * x)
    aff_poly = polynomial((a + 1) * x * y)
    aff_poly_x = polynomial((a + 1) * x)
    @test_throws ErrorException @variable m p Poly(X) unknown_kw = 1
    @test_throws ErrorException @variable m p == 1 Poly(X)
    @test_throws ErrorException @variable m p Poly(X) start = 1
    @test_throws ErrorException @variable m 0 <= p <= 1 Poly(X)
    @test_throws ErrorException @variable m 0 >= p >= 1 Poly(X)
    @test_throws ErrorException @variable m p <= 0 Poly(X)
    @test_throws ErrorException @variable m p >= 1 Poly(X)
    @test_throws ErrorException @variable m p == 1 Poly(X)
    @test_throws ErrorException @variable m p >= 0 Poly(X)
    @test_throws ErrorException @variable(m, p3[2:3] >= 0, Poly(X))
end

function _algebra_element_type(B, mono)
    M = MP.monomial_type(mono)
    C = JuMP.VariableRef
    return SA.AlgebraElement{
        MB.Algebra{MB.SubBasis{B,M,MP.monomial_vector_type(M)},B,M},
        C,
        Vector{C},
    }
end

function _test_variable(
    m,
    p,
    monos,
    binary = false,
    integer = false,
    use_y = true,
    B = MB.Monomial,
)
    @test isa(
        p,
        _algebra_element_type(
            B,
            use_y ? first(monos) : first(variables(first(monos))),
        ),
    )
    @test SA.basis(p).monomials == monomial_vector(monos)
    coeffs = SA.coeffs(p, MB.explicit_basis(p))
    @test all(α -> JuMP.is_binary(α) == binary, coeffs)
    @test all(α -> JuMP.is_integer(α) == integer, coeffs)
end

function test_MonomialBasis(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    X = [1, y, x^2, x, y^2]
    m = Model()
    @variable m p1[1:3] Poly(X)
    PT = _algebra_element_type(MB.Monomial, x * y)
    @test isa(p1, Vector{PT})
    _test_variable(m, p1[1], X)
    @variable(m, p2, Poly(X), integer = true)
    _test_variable(m, p2, X, false, true)
    @variable(m, p3[2:3], Poly(X))
    @test isa(p3, JuMP.Containers.DenseAxisArray{PT,1,Tuple{UnitRange{Int}}})
    _test_variable(m, p3[2], X)
    @variable(m, p4[i = 2:3, j = i:4], Poly(X), binary = true)
    _test_variable(m, p4[2, 3], X, true)

    X = [x^2, y^2]
    @variable m p5[1:3] Poly(X)
    @test isa(p5, Vector{PT})
    _test_variable(m, p5[1], X)
    @variable(m, p6, Poly(MB.SubBasis{MB.Monomial}(X)), integer = true)
    return _test_variable(m, p6, X, false, true)
end

function test_ScaledMonomialBasis(var)
    x = MP.similar_variable(var, Val{:x})
    m = Model()
    @variable(m, p1, Poly(MB.SubBasis{MB.ScaledMonomial}([1, x, x^2])), Int)
    return _test_variable(
        m,
        p1,
        monomial_vector([1, x, x^2]),
        false,
        true,
        false,
        MB.ScaledMonomial,
    )
end

# TODO recover with other basis
#function test_FixedPolynomialBasis(var)
#    x = MP.similar_variable(var, Val{:x})
#    y = MP.similar_variable(var, Val{:y})
#    m = Model()
#    @variable(m, p1, Poly(MB.FixedPolynomialBasis([1 - x^2, x^2 + 2])), Bin)
#    _test_variable(m, p1, monomial_vector([x^2, 1]), true, false, false, false)
#    @variable(m, p2[1:2], Poly(MB.FixedPolynomialBasis([1 - x^2, x^2 + 2])))
#    _test_variable(
#        m,
#        p2[1],
#        monomial_vector([x^2, 1]),
#        false,
#        false,
#        false,
#        false,
#    )
#    # Elements of the basis have type monomial
#    @variable(m, p3[2:3], Poly(MB.FixedPolynomialBasis([x, x^2])))
#    _test_variable(
#        m,
#        p3[2],
#        monomial_vector([x^2, x]),
#        false,
#        false,
#        false,
#        false,
#    )
#    # Elements of the basis have type term
#    @variable(m, p4[1:2], Poly(MB.FixedPolynomialBasis([1, x, x^2])), Int)
#    _test_variable(
#        m,
#        p4[1],
#        monomial_vector([x^2, x, 1]),
#        false,
#        true,
#        false,
#        false,
#    )
#    # Elements of the basis have type variable
#    @variable(
#        m,
#        p5[-1:1],
#        Poly(MB.FixedPolynomialBasis([x, y])),
#        integer = true
#    )
#    _test_variable(m, p5[0], monomial_vector([x, y]), false, true, false)
#    @variable(m, p6[-1:1], Poly(MB.FixedPolynomialBasis([x])), integer = true)
#    return _test_variable(
#        m,
#        p6[0],
#        monomial_vector([x]),
#        false,
#        true,
#        true,
#        false,
#    )
#end

function test_value_function(var)
    x = MP.similar_variable(var, Val{:x})
    y = MP.similar_variable(var, Val{:y})
    m = Model()
    @variable m α
    @variable m β
    p = α * x * y + β * x^2
    JuMP.fix(α, 2)
    JuMP.fix(β, 3)
    expected = 2x * y + 3x^2
    @test_broken JuMP.value(p) == expected
    @test JuMP.value(fix_value, p) == expected
    a = MB.algebra_element(p)
    exp_a = MB.algebra_element(expected)
    @test_broken JuMP.value(a) == exp_a
    @test JuMP.value(fix_value, a) == exp_a
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
TestVariable.runtests(x)
import TypedPolynomials
TypedPolynomials.@polyvar(y)
TestVariable.runtests(y)
