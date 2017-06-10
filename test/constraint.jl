@testset "@constraint macro with polynomials" begin
    m = Model()
    setpolymodule!(m, TestPolyModule)
    @variable m α
    @variable m β
    @polyvar x y
    p = α * x*y + β * x^2
    q = MatPolynomial([α β; β α], [x])
    @test macroexpand(:(@constraint(m, p))).head == :error
    @test macroexpand(:(@constraint(m, begin p >= 0 end))).head == :error
    @test macroexpand(:(@constraint(m, +(p, p, p)))).head == :error
    @test macroexpand(:(@constraint(m, p >= 0, 1))).head == :error
    @test_throws MethodError @constraint(m, p >= 0, unknown_kw=1)
    #@test macroexpand(:(@constraint(m, p >= 0, domain = (@set x >= -1 && x <= 1, domain = y >= -1 && y <= 1)))).head == :error
    @test macroexpand(:(@constraint(m, p + 0, domain = (@set x >= -1 && x <= 1)))).head == :error

    # TODO Once JuMP drops Julia v0.5, this should be move to JuMP and be renamed Base.iszero
    function affexpr_iszero(affexpr)
        affexpr.constant == 0 || return false
        tmp = JuMP.IndexedVector(Float64, m.numCols)
        JuMP.collect_expr!(m, tmp, affexpr)
        tmp.nnz == 0
    end
    _iszero(p) = all(affexpr_iszero, p.a)
    _iszero(p::AbstractArray) = all(_iszero, p)
    function testcon(m, cref, nonnegative, p, ineqs, eqs)
        @test isa(cref, ConstraintRef{Model, PolyJuMP.PolyConstraint})
        c = PolyJuMP.getpolyconstr(m)[cref.idx]
        @test c.nonnegative == nonnegative
        # == between JuMP affine expression is not accurate, e.g. β + α != α + β
        # == 0 is not defined either
        # c.p and p can be matrices
        @test _iszero(c.p - p)
        if isempty(ineqs)
            if isempty(eqs)
                @test isa(c.domain, FullSpace)
            else
                @test isa(c.domain, AlgebraicSet)
                @test c.domain.p == eqs
            end
        else
            @test isa(c.domain, BasicSemialgebraicSet)
            @test c.domain.p == ineqs
            @test c.domain.V.p == eqs
        end
    end

    f(x, y) = @set x + y == 2
    dom = @set x^2 + y^2 == 1 && x^3 + x*y^2 + y >= 1
    testcon(m, @constraint(m, p >= q + 1, domain = @set y >= 1 && dom), true, p - q - 1, [y-1, x^3 + x*y^2 + y - 1], [x^2 + y^2 - 1])
    testcon(m, @constraint(m, p <= q), true, q - p, [], [])
    testcon(m, @constraint(m, p + q >= 0, domain = @set x == y^3), true, p + q, [], [x - y^3])
    testcon(m, @constraint(m, p == q, domain = @set x == 1 && f(x, y)), false, p - q, [], [x - 1, x + y - 2])
    testcon(m, @SDconstraint(m, [p q; q 0] ⪰ [0 0; 0 p]), true, [p q; q -p], [], [])
end

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
    @test macroexpand(:(@polyconstraint(m, p + 0, domain = x >= -1 && x <= 1))).head == :error

    function testcon(m, cref, nonnegative, p, ineqs, eqs)
        @test isa(cref, ConstraintRef{Model, PolyJuMP.PolyConstraint})
        c = PolyJuMP.getpolyconstr(m)[cref.idx]
        @test c.nonnegative == nonnegative
        # == between JuMP affine expression is not accurate, e.g. β + α != α + β
        # == 0 is not defined either
        function affexpr_iszero(affexpr)
            affexpr.constant == 0 || return false
            tmp = JuMP.IndexedVector(Float64, m.numCols)
            JuMP.collect_expr!(m, tmp, affexpr)
            tmp.nnz == 0
        end
        @test all(affexpr_iszero, (c.p - p).a)
        if isempty(ineqs)
            if isempty(eqs)
                @test isa(c.domain, FullSpace)
            else
                @test isa(c.domain, AlgebraicSet)
                @test c.domain.p == eqs
            end
        else
            @test isa(c.domain, BasicSemialgebraicSet)
            @test c.domain.p == ineqs
            @test c.domain.V.p == eqs
        end
    end

    f(x, y) = @set x + y == 2
    dom = @set x^2 + y^2 == 1 && x^3 + x*y^2 + y >= 1
    testcon(m, @polyconstraint(m, p ⪰ q + 1, domain = y >= 1 && dom), true, p - q - 1, [y-1, x^3 + x*y^2 + y - 1], [x^2 + y^2 - 1])
    testcon(m, @polyconstraint(m, p ⪯ q), true, q - p, [], [])
    testcon(m, @polyconstraint(m, p + q >= 0, domain = x == y^3), true, p + q, [], [x - y^3])
    testcon(m, @polyconstraint(m, p == q, domain = x == 1 && f(x, y)), false, p - q, [], [x - 1, x + y - 2])
end
