using Gmsh: Gmsh, gmsh

include("files.jl")

startGmsh()
ShowInfo(1)
rve = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 
# inc1=Cylinder([-0.1, 0.5, 0.5],[2.0,0.0,0.0], 0.1, 2)
# inc2=Sphere([0.5, 0.5, 0.5], 0.25, 2)
# inc3=Cylinder([0.5, -0.1, 0.5],[0.0,2.0,0.0], 0.1, 2)
# inc4=Cylinder([0.5, 0.5,-0.1],[0.0,0.0,2.0], 0.1, 2)

# # inclusions=(Union([inc1,inc2,inc3, inc4],2),)
# inclusions=(Intersection(inc1,inc2,2),)

geo1=Cylinder([-0.1, 0.5, 0.5],[2.0,0.0,0.0], 0.1, 2)
geo2=Cylinder([0.5, -0.1, 0.5],[0.0,2.0,0.0], 0.1, 2)
geo3=Cylinder([0.5, 0.5,-0.1],[0.0,0.0,2.0], 0.1, 2)
geo4=Union([geo1,geo2,geo3],2)

geo5=Sphere([0.5, 0.5, 0.5], 0.25, 2)
L=0.38
geo6=Box([L, L, L], [0.5-L/2, 0.5-L/2, 0.5-L/2], 2)

geo7=Intersection(geo5,geo6,2)
geo8=Cut(geo7,geo4,2)

inclusions=(geo8,)



  


    gmsh.model.add("test")
    xmin, ymin, zmin, xmax, ymax, zmax = _getBoundingBox(rve)
    _addBoundingBox!(gmsh.model, rve)

    dim3::Int32 = 3
    numinc = length(inclusions)
    _tags = Vector{Vector{Tuple{Int32,Int32}}}()
    push!(_tags, [(3, 1)])

    _numinc = zeros(Int32, numinc)
 
 
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









