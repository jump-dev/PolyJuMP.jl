@testset "@variable macro with Poly" begin
    @polyvar x y
    X = [1, y, x^2, x, y^2]

    m = Model()
    @test_throws ErrorException @variable m p Poly(X) unknown_kw=1
    @test_throws ErrorException @variable m p == 1 Poly(X)
    @test_throws ErrorException @variable m p Poly(X) start=1
    @test_throws ErrorException @variable m 0 <= p <= 1 Poly(X)
    @test_throws ErrorException @variable m 0 >= p >= 1 Poly(X)
    @test_throws ErrorException @variable m p <= 0 Poly(X)
    @test_throws ErrorException @variable m p >= 1 Poly(X)
    @test_throws ErrorException @variable m p == 1 Poly(X)
    @test_throws ErrorException @variable m p >= 0 Poly(X)
    @test_throws ErrorException @variable(m, p3[2:3] >= 0, Poly(X))

    function testvar(m, p, x, binary=false, integer=false, vars=true)
        @test isa(p, vars ? DynamicPolynomials.Polynomial{true,JuMP.Variable} : DynamicPolynomials.Polynomial{true,JuMP.AffExpr})
        @test p.x == x
        if vars
            @test all(α -> haskey(m.variabletozeroone, JuMP.index(α)) == binary, coefficients(p))
            @test all(α -> haskey(m.variabletointegrality, JuMP.index(α)) == integer, coefficients(p))
        else
            @test all(α -> haskey(m.variabletozeroone, JuMP.index(first(keys(α.terms)))) == binary, coefficients(p))
            @test all(α -> haskey(m.variabletointegrality, JuMP.index(first(keys(α.terms)))) == integer, coefficients(p))
        end
    end

    @testset "MonomialBasis" begin
        m = Model()
        @variable m p1[1:3] Poly(X)
        @test isa(p1, Vector{DynamicPolynomials.Polynomial{true,JuMP.Variable}})
        testvar(m, p1[1], X)
        @variable(m, p2, Poly(X), integer=true)
        testvar(m, p2, X, false, true)
        @variable(m, p3[2:3], Poly(X))
        @test isa(p3, JuMP.JuMPArray{DynamicPolynomials.Polynomial{true,JuMP.Variable},1,Tuple{UnitRange{Int}}})
        testvar(m, p3[2], X)
        @variable(m, p4[i=2:3,j=i:4], Poly(X), binary=true)
        testvar(m, p4[2,3], X, true)

        X = [x^2, y^2]
        @variable m p5[1:3] Poly(X)
        @test isa(p5, Vector{DynamicPolynomials.Polynomial{true,JuMP.Variable}})
        testvar(m, p5[1], X)
        @variable(m, p6, Poly(X), integer=true)
        testvar(m, p6, X, false, true)
    end

    @testset "FixedPolynomialBasis" begin
        m = Model()
        @variable(m, p1, Poly(FixedPolynomialBasis([1 - x^2, x^2 + 2])), Bin)
        testvar(m, p1, monovec([x^2, 1]), true, false, false)
        @variable(m, p2[1:2], Poly(FixedPolynomialBasis([1 - x^2, x^2 + 2])))
        testvar(m, p2[1], monovec([x^2, 1]), false, false, false)
        # Elements of the basis have type monomial
        @variable(m, p3[2:3], Poly(FixedPolynomialBasis([x, x^2])))
        testvar(m, p3[2], monovec([x^2, x]), false, false, false)
        # Elements of the basis have type term
        @variable(m, p4[1:2], Poly(FixedPolynomialBasis([1, x, x^2])), Int)
        testvar(m, p4[1], monovec([x^2, x, 1]), true, false, false)
        # Elements of the basis have type variable
        @variable(m, p5[-1:1], Poly(FixedPolynomialBasis([x, y])), integer=true)
        testvar(m, p5[0], monovec([x, y]), false, true, false)
        @variable(m, p6[-1:1], Poly(FixedPolynomialBasis([x])), integer=true)
        testvar(m, p6[0], monovec([x]), false, true, false)
    end
end

@testset "getvalue function" begin
    m = Model()
    @variable m α
    @variable m β
    @polyvar x y
    p = α * x*y + β * x^2
    JuMP.fix(α, 2)
    JuMP.fix(β, 3)
    @test getvalue(p) == 2x*y + 3x^2
end
