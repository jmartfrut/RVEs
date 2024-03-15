using MimosaRVE
using Test

@testset "MimosaRVEs.jl" begin

myRVE=RVE([2,2,2],[1.0,1.0,1.0],[1,1,1],[0.0,0.0,0.0])
myInclusion=Sphere([1,1,1],0.25)
model=createGmshModel(myRVE,myInclusion,"test")
createMesh!(model,myRVE, 0.1)
output_file = joinpath(dirname(@__FILE__), "rve.msh")
saveMesh(model, output_file)
# visualizeMesh()
closeGmshModel()

end