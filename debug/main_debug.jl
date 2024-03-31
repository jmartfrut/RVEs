
include("../src/RVEs.jl")
using Main.RVEs: RVEs

RVEs.startGmsh()
RVEs.ShowInfo(0)
myRVE = RVEs.RVE( [1.0, 1.0, 1.0], [0,0,0], [0.0, 0.0, 0.0], 0.1)
myInclusions = RVEs.Sphere([0.5, 0.5, 0.5], 0.2, 2)
model, tagvol=RVEs.createGmshModel(myRVE,myInclusions,"Test1")
RVEs.createMesh!(model,myRVE,myInclusions,tagvol)
RVEs.visualizeMesh()
RVEs.stopGmsh()
 


rve = RVE( [1.0, 1.0, 1.0], [0,0,0], [0.0, 0.0, 0.0], 0.1)
dx=0.05
a=[0.0,0.0,0.0]+[-dx,-dx,-dx]
b=[1.0,1.0,1.0]+[dx,dx,dx]
d=sqrt(sum((b-a).^2))
geo1=Cylinder(a,2*(b-a)/d, 0.05, 2)
a=[0.0,1.0,0.0]+[-dx,dx,-dx]
b=[1.0,0.0,1.0]+[dx,-dx,dx]
d=sqrt(sum((b-a).^2))
geo2=Cylinder(a,2*(b-a)/d, 0.05, 2)
 a=[1.0,0.0,0.0]+[dx,-dx,-dx]
b=[0.0,1.0,1.0]+[-dx,dx,dx]
d=sqrt(sum((b-a).^2))
geo3=Cylinder(a,2*(b-a)/d, 0.05, 2)
a=[1.0,1.0,0.0]+[dx,dx,-dx]
b=[0.0,0.0,1.0]+[-dx,-dx,dx]
d=sqrt(sum((b-a).^2))
geo4=Cylinder(a,2*(b-a)/d, 0.05, 2)
geo5=Sphere([0.5, 0.5, 0.5], 0.2, 2)

geo6=Sphere([0.0, 0.0, 0.0], 0.2, 2)
geo7=Sphere([1.0, 0.0, 0.0], 0.2, 2)
geo8=Sphere([1.0, 1.0, 0.0], 0.2, 2)
geo9=Sphere([0.0,1.0,0.0], 0.2, 2)

geo10=Sphere([0.0, 0.0, 1.0], 0.2, 2)
geo11=Sphere([1.0, 0.0, 1.0], 0.2, 2)
geo12=Sphere([1.0, 1.0, 1.0], 0.2, 2)
geo13=Sphere([0.0,1.0,1.0], 0.2, 2)

geo14=Fuse([geo1, geo2, geo3, geo4, geo5, geo6,geo7,geo8,geo9, geo10,geo11,geo12,geo13],2)

inclusions=(geo14,)

  


    gmsh.model.add("test")
    xmin, ymin, zmin, xmax, ymax, zmax = _getBoundingBox(rve)
    _addBoundingBox!(gmsh.model, rve)

    dim3::Int32 = 3
    numinc = length(inclusions)
    _tags = Vector{Vector{Tuple{Int32,Int32}}}()
    push!(_tags, [(3, 1)])

    _numinc = zeros(Int32, numinc)

    # tag = _addInclusion!(gmsh.model, inclusions[1])

    # gmsh.model.occ.synchronize()
    # visualizeMesh()

 
    for i in 1:numinc
        tag = _addInclusion!(gmsh.model, inclusions[i])
        inc_boundingbox = gmsh.model.occ.getBoundingBox(dim3, tag)
        isinboundary, boundarytype = _isinboundary(rve, inc_boundingbox)
        if isinboundary
            println("Inclusion $i is in the boundary!")
            __tags = _addPeriodicInclusions!(gmsh.model, rve, [(dim3, tag)], boundarytype)
            _numinc[i] = length(__tags)
            push!(_tags, __tags)
        else
            push!(_tags, [(dim3, tag)])
            _numinc[i] = 1
        end
    end

    out_ = []
     _, out__ = gmsh.model.occ.fragment(_tags[1], vcat(_tags[2]...))
     out_matrix, out_inc= _findMatrixVolume(out__)   
     push!(out_, out_inc)

      for i in 2:numinc
          _, out__ = gmsh.model.occ.fragment(out_matrix, vcat(_tags[i+1]...))
        out_matrix, out_inc= _findMatrixVolume(out__)
        @show    out_inc
        if length(out__[1])>1
          push!(out_, out_inc)
        end
      end
     push!(out_, out_matrix)
  
     gmsh.model.occ.synchronize()
     removol = gmsh.model.getEntities(3)

    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    eps = 1e-3
    vin = gmsh.model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 3)
    for v in vin
        deleteat!(removol, findall(x -> x == v, removol))
    end
    gmsh.model.removeEntities(removol, true)

     _addPhysicalGroups!(gmsh.model, out_, rve)

     createMesh!(gmsh.model,rve,inclusions,out_)

     visualizeMesh()









