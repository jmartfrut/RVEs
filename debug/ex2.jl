using Gmsh: Gmsh, gmsh

include("files.jl")

startGmsh()
ShowInfo(1)
rve = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 numfib=[8,8] # y, z
 step = rve.size[2:3]/(numfib[1]-1)
 gamma = 0.1*step
 radius = minimum((step/2.0)-gamma)

 myInclusions=[]
for iz in 1:numfib[2]
  for iy in 1:numfib[1]
    z = (iz-1)*step[2]
    y = (iy-1)*step[1]
    push!(myInclusions, Cylinder([-0.1, y, z],[2.0,0.0,0.0], radius, 2))
  end
end

inclusions=tuple(myInclusions...)
# rve = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
# inclusions = (Sphere([0.0, 0.5, 0.5], 0.25, 2),
# Sphere([0.5, 0.0, 0.5], 0.25, 2),
# Sphere([0.5, 0.5, 0.0], 0.25, 2),
# Sphere([1.0, 1.0, 1.0], 0.25, 2))



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









