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

    function testvar(p::MultivariatePolynomials.AbstractPolynomial, nonnegative, monotype, category=:Cont)
        @test !nonnegative
    end
    function testvar(p::TestPolyModule.TestPolyVar{MS}, nonnegative, monotype, category=:Cont) where MS
        @test MS == monotype
        @test p.x == monovec(X)
        @test p.x[1] == x^2
        @test p.x[2] == y^2
        @test p.x[3] == x
        @test p.x[4] == y
        @test p.x[5] == 1
        @test p.category == category
    end

    @variable m p1[1:3] Poly{true}(X)
    @test isa(p1, Vector{TestPolyModule.TestPolyVar{:Default}})
    testvar(p1[1], true, :Default)
    @variable(m, p2, Poly{false, :Classic}(X), category=:Int)
    testvar(p2, false, :Classic, :Int)
    @variable(m, p3, Poly{false, :Gram}(X))
    testvar(p3, false, :Gram)
    @variable m p4 >= 0 Poly{true}(X)
    testvar(p4, true, :Default)
    @variable(m, p5 >= 0, Poly{true, :Classic}(X))
    testvar(p5, true, :Classic)
    @variable(m, p6[2:3] >= 0, Poly{true, :Gram}(X))
    @test isa(p6, JuMP.JuMPArray{TestPolyModule.TestPolyVar{:Gram},1,Tuple{UnitRange{Int}}})
    testvar(p6[2], true, :Gram)
    @variable(m, p7[i=2:3,j=i:4], Poly{true}(X), category=:Bin)
    testvar(p7[2,3], true, :Default, :Bin)
end

@testset "@variable macro with Poly: Default methods" begin
    m = Model()
    @polyvar x y
    X = [x^2, y^2]

    function testvar(p, nonnegative, monotype, x, category=:Cont)
        @test isa(p, DynamicPolynomials.Polynomial{true,JuMP.Variable})
        @test p.x == x
    end

    @variable m p1[1:3] Poly(X)
    @test isa(p1, Vector{DynamicPolynomials.Polynomial{true,JuMP.Variable}})
    testvar(p1[1], false, :Default, X)
    @variable(m, p2, Poly{false, :Classic}(X), category=:Int)
    testvar(p2, false, :Classic, X, :Int)
    @variable(m, p3, Poly{false, :Gram}(X))
    testvar(p3, false, :Gram, [x^4,x^2*y^2,y^4])
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
