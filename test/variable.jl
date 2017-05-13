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
    @test macroexpand(:(@polyvariable(m, p, gramonomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, grammonomials=X, monomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, monomials=X, grammonomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, X, monomials=X))).head == :error
    @test macroexpand(:(@polyvariable(m, p, X, X))).head == :error

    function testvar(p, nonnegative, monotype, x)
        @test isa(p, TestPolyModule.TestPoly)
        @test p.nonnegative == nonnegative
        @test p.monotype == monotype
        @test p.x == x
    end

    @polyvariable m p1 X
    testvar(p1, false, :Default, X)
    @polyvariable(m, p2, monomials=X)
    testvar(p2, false, :Classic, X)
    @polyvariable(m, p3, grammonomials=X)
    testvar(p3, false, :Gram, X)
    @polyvariable m p4 >= 0 X
    testvar(p4, true, :Default, X)
    @polyvariable(m, p5 >= 0, monomials=X)
    testvar(p5, true, :Classic, X)
    @polyvariable(m, p6 >= 0, grammonomials=X)
    testvar(p6, true, :Gram, X)
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
