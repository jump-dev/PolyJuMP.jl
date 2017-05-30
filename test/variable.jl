@testset "@variable macro with Poly" begin
    m = Model()
    setpolymodule!(m, TestPolyModule)
    @polyvar x y
    X = [x^2, y^2]

    #@test_throws ErrorException @variable m p Poly(X) unknown_kw=1 # Fails on v0.5 on Travis
    @test_throws ErrorException @variable m p == 1 Poly(X)
    @test_throws ErrorException @variable m 0 <= p <= 1 Poly(X)
    @test_throws ErrorException @variable m 0 >= p >= 1 Poly(X)
    @test_throws ErrorException @variable m p <= 0 Poly(X)
    @test_throws ErrorException @variable m p >= 1 Poly(X)
    @test_throws ErrorException @variable m p == 1 Poly(X)
    @test_throws ErrorException @variable m p >= 0 Poly(X)

    function testvar(p, nonnegative, monotype, x, category=:Cont)
        @test isa(p, TestPolyModule.TestPoly{nonnegative})
        @test p.monotype == monotype
        @test p.x == x
        @test p.category == category
    end

    @variable m p1[1:3] Poly(X)
    @test isa(p1, Vector{TestPolyModule.TestPoly{false}})
    testvar(p1[1], false, :Default, X)
    @variable(m, p2, Poly{false, :Classic}(X), category=:Int)
    testvar(p2, false, :Classic, X, :Int)
    @variable(m, p3, Poly{false, :Gram}(X))
    testvar(p3, false, :Gram, X)
    @variable m p4 >= 0 Poly{true}(X)
    testvar(p4, true, :Default, X)
    @variable(m, p5 >= 0, Poly{true, :Classic}(X))
    testvar(p5, true, :Classic, X)
    @variable(m, p6[2:3] >= 0, Poly{true, :Gram}(X))
    @test isa(p6, JuMP.JuMPArray{TestPolyModule.TestPoly{true},1,Tuple{UnitRange{Int64}}})
    testvar(p6[2], true, :Gram, X)
    @variable(m, p7[i=2:3,j=i:4], Poly{true}(X), category=:Bin)
    testvar(p7[2,3], true, :Default, X, :Bin)
end

@testset "getvalue function" begin
    m = Model()
    @variable m α
    @variable m β
    @polyvar x y
    p = α * x*y + β * x^2
    q = MatPolynomial([α β; β α], [x, y])
    JuMP.fix(α, 2)
    JuMP.fix(β, 3)
    @test getvalue(p) == 2x*y + 3x^2
    # Explicit polynomial conversion is needed only if MultivariatePolynomials < v0.0.2
    @test Polynomial(getvalue(q)) == 2x^2 + 2y^2 + 6x*y
end

@testset "@polyvariable macro" begin
    m = Model()
    setpolymodule!(m, TestPolyModule)
    @polyvar x y
    X = [x^2, y^2]

    @test macroexpand(:(@polyvariable p >= 0 X)).head == :error
    @test macroexpand(:(@polyvariable m 0 <= p <= 1 X)).head == :error
    @test macroexpand(:(@polyvariable m 0 >= p >= 1 X)).head == :error
    @test macroexpand(:(@polyvariable m p <= 0 X)).head == :error
    @test macroexpand(:(@polyvariable m p >= 1 X)).head == :error
    @test macroexpand(:(@polyvariable m p == 1 X)).head == :error
    @test macroexpand(:(@polyvariable m p)).head == :error
    @test macroexpand(:(@polyvariable m p[1:2] X)).head == :error
    @test macroexpand(:(@polyvariable(m, p, gramonomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, grammonomials=X, monomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, monomials=X, grammonomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, X, monomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, X, X))).head == :error

    function testvar(p, nonnegative, monotype, x, category=:Cont)
        @test isa(p, TestPolyModule.TestPoly{nonnegative})
        @test p.monotype == monotype
        @test p.x == x
        @test p.category == category
    end

    @polyvariable m p1 X
    testvar(p1, false, :Default, X)
    @polyvariable(m, p2, monomials=X, category=:Int)
    testvar(p2, false, :Classic, X, :Int)
    @polyvariable(m, p3, grammonomials=X)
    testvar(p3, false, :Gram, X)
    @polyvariable m p4 >= 0 X
    testvar(p4, true, :Default, X)
    @polyvariable(m, p5 >= 0, monomials=X)
    testvar(p5, true, :Classic, X)
    @polyvariable(m, p6 >= 0, grammonomials=X)
    testvar(p6, true, :Gram, X)
    @polyvariable(m, p7 >= 0, X, category=:Bin)
    testvar(p7, true, :Default, X, :Bin)
end
