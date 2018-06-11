@testset "@variable macro with Poly" begin
    m = Model()
    setpolymodule!(m, TestPolyModule)
    @polyvar x y
    X = [1, y, x^2, x, y^2]

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

    function testvar(p, x, category=:Cont)
        @test isa(p, DynamicPolynomials.Polynomial{true,JuMP.Variable})
        @test p.x == x
        @test all(α -> m.colCat[α.col] == category, coefficients(p))
    end

    @variable m p1[1:3] Poly(X)
    @test isa(p1, Vector{DynamicPolynomials.Polynomial{true,JuMP.Variable}})
    testvar(p1[1], X)
    @variable(m, p2, Poly(X), category=:Int)
    testvar(p2, X, :Int)
    @variable(m, p3[2:3], Poly(X))
    @test isa(p3, JuMP.JuMPArray{DynamicPolynomials.Polynomial{true,JuMP.Variable},1,Tuple{UnitRange{Int}}})
    testvar(p3[2], X)
    @variable(m, p4[i=2:3,j=i:4], Poly(X), category=:Bin)
    testvar(p4[2,3], X, :Bin)

    X = [x^2, y^2]
    @variable m p5[1:3] Poly(X)
    @test isa(p5, Vector{DynamicPolynomials.Polynomial{true,JuMP.Variable}})
    testvar(p5[1], X)
    @variable(m, p6, Poly(X), category=:Int)
    testvar(p6, X, :Int)

    @variable(m, p7, Poly(FixedPolynomialBasis([1 - x^2, x^2 + 2])), category=:Bin)
    testvar(p7, monovec([x^2, 1]), :Bin)
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
