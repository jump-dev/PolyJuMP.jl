@testset "PolyModule" begin
    m = Model()
    # Triggers the creation of polydata
    @test PolyJuMP.getpolydata(m).nonnegpoly_default === nothing
    @test_throws ErrorException PolyJuMP.getdefault(m, PolyJuMP.NonNegPoly)
    @test PolyJuMP.getpolydata(m).posdefpolymatrix_default === nothing
    @test_throws ErrorException PolyJuMP.getdefault(m, PolyJuMP.PosDefPolyMatrix)
    setpolymodule!(m, TestPolyModule)
    @test PolyJuMP.getdefault(m, PolyJuMP.NonNegPoly) == TestPolyModule.TestNonNeg
    @test PolyJuMP.getdefault(m, PolyJuMP.PosDefPolyMatrix) == TestPolyModule.TestPosDefMatrix
    PolyJuMP.setdefault!(m, PolyJuMP.PosDefPolyMatrix, TestPolyModule.TestNonNeg)
    @test PolyJuMP.getdefault(m, PolyJuMP.PosDefPolyMatrix) == TestPolyModule.TestNonNeg
end
