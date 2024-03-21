using MimosaRVE


using Test

# @testset "MimosaRVEs.jl" begin



myRVE = RVE([2, 2, 2], [1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
myInclusions = (Sphere([1.0, 1.0, 0.5], 0.25, 0.01), Sphere([0.5, 0.5, 1.0], 0.25, 0.05), Sphere([0.5, 0.5, 0.5], 0.2, 0.03))
model, tagvol=createGmshModel(myRVE,myInclusions)
createMesh!(model,myRVE,myInclusions,tagvol)
output_file = joinpath(dirname(@__FILE__), "rve.msh")
saveMesh(model, output_file)
visualizeMesh()
closeGmshModel()
 