
function _addBoundingBox!(model::Module, rve::RVE)
    model.occ.addBox(rve.origin[1], rve.origin[2], rve.origin[3], rve.size[1], rve.size[2], rve.size[3], 1)
    model.occ.synchronize()
end

function _isinboundary(rve::RVE, b::NTuple{6,Float64})
    xmin, ymin, zmin, xmax, ymax, zmax = _getBoundingBox(rve)
    out = false
    location = [0, 0, 0, 0, 0, 0]
    if b[1] <= xmin
        out = true
        location[1] = 1
    end
    if b[2] <= ymin
        out = true
        location[2] = 1
    end
    if b[3] <= zmin
        out = true
        location[3] = 1
    end
    if b[4] >= xmax
        out = true
        location[4] = 1
    end
    if b[5] >= ymax
        out = true
        location[5] = 1
    end
    if b[6] >= zmax
        out = true
        location[6] = 1
    end

    return out, location
end


function _getBoundingBox(rve::RVE)
    xmin = rve.origin[1]
    ymin = rve.origin[2]
    zmin = rve.origin[3]
    xmax = rve.origin[1] + rve.size[1]
    ymax = rve.origin[2] + rve.size[2]
    zmax = rve.origin[3] + rve.size[3]
    return [xmin, ymin, zmin, xmax, ymax, zmax]
end
 

function _addPhysicalGroups!(model::Module, out_::Vector{Any}, rve::RVE)
    xmin, ymin, zmin, xmax, ymax, zmax = _getBoundingBox(rve)
    eps = 1e-3

    allinc = vcat(out_[1:end-1]...)
    model.addPhysicalGroup(3, [out_[end][2]], 1, "Phase0")
    numinc = length(out_) - 1
    for i in 1:numinc
        model.addPhysicalGroup(3, map(x -> x[2], vcat(out_[i]...)), i + 2, "Phase$i")
    end
    model.addPhysicalGroup(3, map(x -> x[2], allinc), 2, "Inclusions")

    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmin + eps, ymax + eps, zmax + eps, 0)
    ent1 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmin + eps, ymax + eps, zmax + eps, 1)
    ent2 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmin + eps, ymax + eps, zmax + eps, 2)
    model.addPhysicalGroup(0, map(x -> x[2], ent0), 1, "Pointsxmin")
    model.addPhysicalGroup(1, map(x -> x[2], ent1), 1, "Linesxmin")
    model.addPhysicalGroup(2, map(x -> x[2], ent2), 1, "Surfxmin")

    ent0 = model.getEntitiesInBoundingBox(xmax - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 0)
    ent1 = model.getEntitiesInBoundingBox(xmax - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 1)
    ent2 = model.getEntitiesInBoundingBox(xmax - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 2)
    model.addPhysicalGroup(0, map(x -> x[2], ent0), 2, "Pointsxmax")
    model.addPhysicalGroup(1, map(x -> x[2], ent1), 2, "Linesxmax")
    model.addPhysicalGroup(2, map(x -> x[2], ent2), 2, "Surfxmax")

    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymin + eps, zmax + eps, 0)
    ent1 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymin + eps, zmax + eps, 1)
    ent2 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymin + eps, zmax + eps, 2)
    model.addPhysicalGroup(0, map(x -> x[2], ent0), 3, "Pointsymin")
    model.addPhysicalGroup(1, map(x -> x[2], ent1), 3, "Linesymin")
    model.addPhysicalGroup(2, map(x -> x[2], ent2), 3, "Surfymin")

    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymax - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 0)
    ent1 = model.getEntitiesInBoundingBox(xmin - eps, ymax - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 1)
    ent2 = model.getEntitiesInBoundingBox(xmin - eps, ymax - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 2)
    model.addPhysicalGroup(0, map(x -> x[2], ent0), 4, "Pointsymax")
    model.addPhysicalGroup(1, map(x -> x[2], ent1), 4, "Linesymax")
    model.addPhysicalGroup(2, map(x -> x[2], ent2), 4, "Surfymax")


    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmin + eps, 0)
    ent1 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmin + eps, 1)
    ent2 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmin + eps, 2)
    model.addPhysicalGroup(0, map(x -> x[2], ent0), 5, "Pointszmin")
    model.addPhysicalGroup(1, map(x -> x[2], ent1), 5, "Lineszmin")
    model.addPhysicalGroup(2, map(x -> x[2], ent2), 5, "Surfzmin")


    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmax - eps, xmax + eps, ymax + eps, zmax + eps, 0)
    ent1 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmax - eps, xmax + eps, ymax + eps, zmax + eps, 1)
    ent2 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmax - eps, xmax + eps, ymax + eps, zmax + eps, 2)
    model.addPhysicalGroup(0, map(x -> x[2], ent0), 6, "Pointszmax")
    model.addPhysicalGroup(1, map(x -> x[2], ent1), 6, "Lineszmax")
    model.addPhysicalGroup(2, map(x -> x[2], ent2), 6, "Surfzmax")

    corners = []
    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmin + eps, ymin + eps, zmin + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmax - eps, ymin - eps, zmin - eps, xmax + eps, ymin + eps, zmin + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymax - eps, zmin - eps, xmin + eps, ymax + eps, zmin + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmax - eps, ymax - eps, zmin - eps, xmax + eps, ymax + eps, zmin + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmax - eps, xmin + eps, ymin + eps, zmax + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmax - eps, ymin - eps, zmax - eps, xmax + eps, ymin + eps, zmax + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmin - eps, ymax - eps, zmax - eps, xmin + eps, ymax + eps, zmax + eps, 0)
    push!(corners, map(x -> x[2], ent0))
    ent0 = model.getEntitiesInBoundingBox(xmax - eps, ymax - eps, zmax - eps, xmax + eps, ymax + eps, zmax + eps, 0)
    push!(corners, map(x -> x[2], ent0))

    model.addPhysicalGroup(0, vcat(corners...), 7, "Corners")

end


function _findMatrixVolume(out_::Vector{Vector{Tuple{Int32, Int32}}})
    out_raw=out_[1]
    out_raw_=vcat(out_[2:end]...)
    out_matrix=[]
    out_inclusion=[]
    for j in  out_raw
        if  !(j in out_raw_)
            out_matrix=j
        else
            push!(out_inclusion, j)
        end
    end
    return out_matrix, out_inclusion
end


function createGmshModel(rve::RVE, inclusion::Inclusion ,name::String)
    return createGmshModel(rve, (inclusion, ), name)
end


function createGmshModel(rve::RVE, inclusions::Tuple, name::String)

    gmsh.model.add(name)
    xmin, ymin, zmin, xmax, ymax, zmax = _getBoundingBox(rve)
    _addBoundingBox!(gmsh.model, rve)
    dim3::Int32 = 3
    numinc = length(inclusions)
    _tags = Vector{Vector{Tuple{Int32,Int32}}}()
    push!(_tags, [(3, 1)])

    _numinc = zeros(Int32, numinc)
    for i in 1:numinc
        tag = inclusions[i](gmsh.model)
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

    return gmsh.model, out_
end
 


function _addPeriodicInclusions!(model::Module, rve::RVE, _tags::Vector{Tuple{Int32,Int32}}, location::Vector{Int64})


    if rve.periodicityFlags[1] == 1 && (location[1] == 1 || location[4] == 1) && !(location[1]==1 && location[4]==1)
        newtags_x = model.occ.copy(_tags)

        if location[1] == 1
            gmsh.model.occ.translate(newtags_x, rve.size[1], 0.0, 0.0)
        elseif location[4] == 1
            gmsh.model.occ.translate(newtags_x, -rve.size[1], 0.0, 0.0)
        end
        _tags = vcat(_tags, newtags_x)
    end

    if rve.periodicityFlags[2] == 1 && (location[2] == 1 || location[5] == 1)&& !(location[2]==1 && location[5]==1)
        newtags_y = model.occ.copy(_tags)

        if location[2] == 1
            gmsh.model.occ.translate(newtags_y, 0.0, rve.size[2], 0.0)
        elseif location[5] == 1
            gmsh.model.occ.translate(newtags_y, 0.0, -rve.size[2], 0.0)
        end
        _tags = vcat(_tags, newtags_y)
    end

    if rve.periodicityFlags[3] == 1 && (location[3] == 1 || location[6] == 1)&& !(location[3]==1 && location[6]==1)
        newtags_z = model.occ.copy(_tags)

        if location[3] == 1
            gmsh.model.occ.translate(newtags_z, 0.0, 0.0, rve.size[3])
        elseif location[6] == 1
            gmsh.model.occ.translate(newtags_z, 0.0, 0.0, -rve.size[3])
        end
        _tags = vcat(_tags, newtags_z)
    end
    return _tags
end
