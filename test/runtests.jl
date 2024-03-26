using RVEs
using Test

@testset "Test1: Simple Cubic nonperiodic" begin
startGmsh()
ShowInfo(0)
myRVE = RVE( [1.0, 1.0, 1.0], [0,0,0], [0.0, 0.0, 0.0], 0.1)
myInclusions = (Sphere([0.5, 0.5, 0.5], 0.2, 2),)
model, tagvol=createGmshModel(myRVE,myInclusions,"Test1")
createMesh!(model,myRVE,myInclusions,tagvol)
_, nodeCoords, _ = model.mesh.getNodes()
output_file = joinpath(dirname(@__FILE__), "Test1.msh")
saveMesh(output_file)
@test div(length(nodeCoords),3)==1827
@test nodeCoords[1276]==0.5104247300511138
stopGmsh()
end


@testset "Test2: Simple Cubic periodic" begin
  startGmsh()
  ShowInfo(0)
  myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
  myInclusions = (Sphere([0.5, 0.5, 0.5], 0.2, 2),)
  model, tagvol=createGmshModel(myRVE,myInclusions,"Test2")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
  @test div(length(nodeCoords),3)==1835
  @test nodeCoords[2475]==0.26005762047093306
  output_file = joinpath(dirname(@__FILE__), "Test2.msh")
  saveMesh(output_file)
  stopGmsh()
end

  
@testset "Test3: Sphere+Ellipsoid" begin
  startGmsh()
  ShowInfo(0)
 myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 myInclusions = (Sphere([1.0, 1.0, 0.5], 0.1, 2), 
                 Ellipsoid([0.5, 0.5, 0.5], [0.25, 0.1, 0.1], [pi/4, pi/4], 2))
  model, tagvol=createGmshModel(myRVE,myInclusions,"Test3")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
  @test div(length(nodeCoords),3)==3279
  @test nodeCoords[1276]==0.374607510011319
  output_file = joinpath(dirname(@__FILE__), "Test3.msh")
  saveMesh(output_file)
  stopGmsh()
  end


@testset "Test4: Body Centered Cubic" begin
  startGmsh()
  ShowInfo(0)
 myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 myInclusions = (Sphere([0.5, 0.5, 0.5], 0.25, 2),
 Sphere([1.0, 1.0, 1.0], 0.25, 2))
  model, tagvol=createGmshModel(myRVE,myInclusions,"Test4")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
  @test div(length(nodeCoords),3)==3045
  @test nodeCoords[1276]==0.6595175874937982
  output_file = joinpath(dirname(@__FILE__), "Test4.msh")
  saveMesh(output_file)
  stopGmsh()
  end



@testset "Test5: Face Centered Cubic" begin
  startGmsh()
  ShowInfo(0)
 myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 myInclusions = (Sphere([0.0, 0.5, 0.5], 0.25, 2),
 Sphere([0.5, 0.0, 0.5], 0.25, 2),
 Sphere([0.5, 0.5, 0.0], 0.25, 2),
 Sphere([1.0, 1.0, 1.0], 0.25, 2))
  model, tagvol=createGmshModel(myRVE,myInclusions,"Test5")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
  @test div(length(nodeCoords),3)==4778
  @test nodeCoords[1276]==0.951227419495968
  output_file = joinpath(dirname(@__FILE__), "Test5.msh")
  saveMesh(output_file)
  stopGmsh()
  end



@testset "Test6: Face Centered Fibers" begin
  startGmsh()
  ShowInfo(0)
 myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 myInclusions = (Cylinder([-0.1, 0.5, 0.5],[2.0,0.0,0.0], 0.25, 2),
 Cylinder([-0.1, 1.0, 1.0], [2.0,0.0,0.0], 0.25, 2))
  model, tagvol=createGmshModel(myRVE,myInclusions,"Test6")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
  @test div(length(nodeCoords),3)==4963
  @test nodeCoords[2435]==0.3754847230467798
  output_file = joinpath(dirname(@__FILE__), "Test6.msh")
  saveMesh(output_file)
  stopGmsh()
  end



@testset "Test7: Face Centered matrix of Fibers parametric" begin
  startGmsh()
  ShowInfo(0)
 myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 numfib=[4,4] # y, z
 step = myRVE.size[2:3]/(numfib[1]-1)
 gamma = 0.1*step
 radius = minimum((step/2.0)-gamma)

 myInclusions_=[]
for iz in 1:numfib[2]
  for iy in 1:numfib[1]
    z = (iz-1)*step[2]
    y = (iy-1)*step[1]
    push!(myInclusions_, Cylinder([-0.1, y, z],[2.0,0.0,0.0], radius, 2))
  end
end
myInclusions=tuple(myInclusions_...)
  model, tagvol=createGmshModel(myRVE,tuple(myInclusions...),"Test7")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
  # @test div(length(nodeCoords),3)==4963
  # @test nodeCoords[2435]==0.3754847230467798
  output_file = joinpath(dirname(@__FILE__), "Test7.msh")
  saveMesh(output_file)
  stopGmsh()
  end



@testset "Test8: Inclusion from bolean operations" begin
  startGmsh()
  ShowInfo(0)
 myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 geo1=Cylinder([-0.1, 0.5, 0.5],[2.0,0.0,0.0], 0.1, 2)
 geo2=Cylinder([0.5, -0.1, 0.5],[0.0,2.0,0.0], 0.1, 2)
 geo3=Cylinder([0.5, 0.5,-0.1],[0.0,0.0,2.0], 0.1, 2)
 geo4=Fuse([geo1,geo2,geo3],2)
 geo5=Sphere([0.5, 0.5, 0.5], 0.25, 2)
 L=0.38
 geo6=Box([L, L, L], [0.5-L/2, 0.5-L/2, 0.5-L/2], 2)
 geo7=Intersect(geo5,geo6,2)
 geo8=Cut(geo7,geo4,2)
 myInclusions=(geo8,)
  model, tagvol=createGmshModel(myRVE,myInclusions,"Test8")
  createMesh!(model,myRVE,myInclusions,tagvol)
  _, nodeCoords, _ = model.mesh.getNodes()
   @test div(length(nodeCoords),3)==4686
   @test nodeCoords[2435]==0.3215734631563863
  output_file = joinpath(dirname(@__FILE__), "Test8.msh")
  saveMesh(output_file)
  stopGmsh()
  end



 