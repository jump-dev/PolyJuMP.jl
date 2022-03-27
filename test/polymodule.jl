module TestPolyModule

using Test

using JuMP
using PolyJuMP

include("testpolymodule.jl")

function test_PolyModule()
    m = Model()
    # Triggers the creation of polydata
    @test PolyJuMP.getpolydata(m).nonnegpoly_default === nothing
    @test_throws ErrorException PolyJuMP.getdefault(m, PolyJuMP.NonNegPoly)
    @test PolyJuMP.getpolydata(m).posdefpolymatrix_default === nothing
    @test_throws ErrorException PolyJuMP.getdefault(m, PolyJuMP.PosDefPolyMatrix)
    setpolymodule!(m, DummyPolyModule)
    @test PolyJuMP.getdefault(m, PolyJuMP.NonNegPoly) == DummyPolyModule.DummyNonNeg
    @test PolyJuMP.getdefault(m, PolyJuMP.PosDefPolyMatrix) == DummyPolyModule.DummyPosDefMatrix
    PolyJuMP.setdefault!(m, PolyJuMP.PosDefPolyMatrix, DummyPolyModule.DummyNonNeg)
    @test PolyJuMP.getdefault(m, PolyJuMP.PosDefPolyMatrix) == DummyPolyModule.DummyNonNeg
end

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$name", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
end

end

TestPolyModule.runtests()
