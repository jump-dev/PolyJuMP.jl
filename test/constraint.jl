@testset "@polyconstraint macro" begin
    m = Model()
    setpolymodule!(m, TestPolyModule)
    @variable m α
    @variable m β
    @polyvar x y
    p = α * x*y + β * x^2
    q = MatPolynomial([α β; β α], [x])
    @test macroexpand(:(@polyconstraint(m, p))).head == :error
    @test macroexpand(:(@polyconstraint(m, begin p >= 0 end))).head == :error
    @test macroexpand(:(@polyconstraint(m, +(p, p, p)))).head == :error
    @test macroexpand(:(@polyconstraint(m, p >= 0, 1))).head == :error
    @test macroexpand(:(@polyconstraint(m, p >= 0, unknown_kw=1))).head == :error
    @test macroexpand(:(@polyconstraint(m, p >= 0, domain = x >= -1 && x <= 1, domain = y >= -1 && y <= 1))).head == :error
    @test macroexpand(:(@polyconstraint(m, p + 0, domain = x >= -1 && x <= 1, domain = y >= -1 && y <= 1))).head == :error

    function testcon(m, cref, nonnegative, p, ineqs, eqs)
        c = PolyJuMP.getpolyconstr(m)[cref.idx]
        @test isa(c, TestPolyModule.TestCon)
        @test c.nonnegative == nonnegative
        @test c.p == p
        @show typeof(c.domain)
    end

    #testcon(m, @polyconstraint(m, p ⪰ q + 1, domain = y >= 1 && x^2 + y^2 == 1 && x^3 + x*y^2 + y >= 1), true, p - q - 1, [], [])
    @polyconstraint(m, p ⪯ q)
    @polyconstraint(m, p + q >= 0, domain = x == y^3)
    @polyconstraint(m, p == q, domain = x == 1 && x + y == 2)
end
