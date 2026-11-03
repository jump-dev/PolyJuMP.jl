module TestBridgeable

using Test
using JuMP
using PolyJuMP
import MathOptInterface as MOI

# These sets deliberately list bridges that do not support them.
struct UnsupportedSet{B} <: MOI.AbstractVectorSet end

function PolyJuMP.bridges(::Type{UnsupportedSet{B}}) where {B}
    return Tuple{Type,Type}[(B, Float64)]
end

function PolyJuMP.bridges(
    ::Type{<:MOI.AbstractVectorFunction},
    S::Type{<:UnsupportedSet},
)
    return PolyJuMP.bridges(S)
end

function test_supports()
    VB = MOI.Bridges.Variable.ZerosBridge{Float64}
    @test PolyJuMP._supports(VB, MOI.Zeros)
    @test !PolyJuMP._supports(VB, MOI.Nonnegatives)
    CB = MOI.Bridges.Constraint.ScalarizeBridge{Float64}
    @test PolyJuMP._supports(CB, MOI.Zeros)
    @test !PolyJuMP._supports(CB, MOI.SecondOrderCone)
end

function test_unsupported_bridges()
    c = JuMP.VectorConstraint([1.0], MOI.Zeros(1))
    VS = UnsupportedSet{MOI.Bridges.Variable.ZerosBridge}
    @test_throws AssertionError PolyJuMP.bridgeable(c, VS)
    CS = UnsupportedSet{MOI.Bridges.Constraint.ScalarizeBridge}
    @test_throws AssertionError PolyJuMP.bridgeable(c, CS)
    @test_throws AssertionError PolyJuMP.bridgeable(
        c,
        MOI.VectorOfVariables,
        CS,
    )
end

@testset "Bridge support" begin
    test_supports()
    test_unsupported_bridges()
end

end
