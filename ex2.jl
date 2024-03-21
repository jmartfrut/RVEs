using Gmsh: Gmsh, gmsh

abstract type Inclusion end

struct RVE
    nCells::Vector{Int64}    # number of cells in each direction
    size::Vector{Float64}    # size of the RVE
    periodicityFlags::Vector{Int64}       # periodicity flags
    origin::Vector{Float64}     # set cell origin to [0,0,0]
    meshsize::Float64
end

struct Sphere <: Inclusion
    origin::Vector{Float64} # center of the sphere
    radius::Float64 # radius of the sphere
    meshsize::Float64
end


struct Cylinder <: Inclusion
    origin::Vector{Float64} # center of the sphere
    axis::Vector{Float64}  # axis of cilinder
    radius::Float64 #radius of cylinder
    meshsize::Float64
end


function _addInclusion!(model::Module, inc::Inclusion)
    if inc isa Sphere
        tag = model.occ.addSphere(inc.origin[1], inc.origin[2], inc.origin[3], inc.radius, -1)
    elseif inc isa Cylinder
        tag = model.occ.addCylinder(inc.origin[1], inc.origin[2], inc.origin[3], inc.axis[1], inc.axis[2], inc.axis[3], inc.radius, -1)
    end
    return tag
end

function _addBoundingBox!(model::Module, rve::RVE)
    xmin, ymin, zmin, xmax, ymax, zmax = getBoundingBox(rve::RVE)
    model.occ.addBox(xmin, ymin, zmin, xmax, ymax, zmax, 1)
    model.occ.synchronize()
end

function _isinboundary(rve::RVE, b::NTuple{6,Float64})
    xmin, ymin, zmin, xmax, ymax, zmax = getBoundingBox(rve)
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

function _addPeriodicInclusions!(model::Module, rve::RVE, _tags::Vector{Tuple{Int32,Int32}}, location::Vector{Int64})


    if rve.periodicityFlags[1] == 1 && (location[1] == 1 || location[4] == 1)
        newtags_x = model.occ.copy(_tags)

        if location[1] == 1
            gmsh.model.occ.translate(newtags_x, rve.size[1], 0.0, 0.0)
        elseif location[4] == 1
            gmsh.model.occ.translate(newtags_x, -rve.size[1], 0.0, 0.0)
        end
        _tags = vcat(_tags, newtags_x)
    end

    if rve.periodicityFlags[2] == 1 && (location[2] == 1 || location[5] == 1)
        newtags_y = model.occ.copy(_tags)

        if location[2] == 1
            gmsh.model.occ.translate(newtags_y, 0.0, rve.size[2], 0.0)
        elseif location[5] == 1
            gmsh.model.occ.translate(newtags_y, 0.0, -rve.size[2], 0.0)
        end
        _tags = vcat(_tags, newtags_y)
    end

    if rve.periodicityFlags[3] == 1 && (location[3] == 1 || location[6] == 1)
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

function createGmshModel(rve::RVE, inclusions::Tuple, namedModel::String)
    gmsh.initialize()
    gmsh.model.add(namedModel)
    xmin, ymin, zmin, xmax, ymax, zmax = getBoundingBox(rve)
    _addBoundingBox!(gmsh.model, rve)
    dim3::Int32 = 3
    numinc = length(inclusions)
    _tags::Array{Tuple{Int32,Int32},1} = []
    for i in 1:1
        tag = _addInclusion!(gmsh.model, inclusions[i])
        gmsh.model.occ.synchronize()
        inc_boundingbox = gmsh.model.getBoundingBox(dim3, tag)
        isinboundary = _isinboundary(rve, inc_boundingbox)
        _tags = vcat(_tags, [(dim3, tag)])
        if isinboundary
            println("Inclusion $i is in the boundary!")
            _tags = _addPeriodicInclusions!(gmsh.model, rve, _tags)
        end
    end

    out, _ = gmsh.model.occ.fragment([(3, 1)], _tags)
    gmsh.model.occ.synchronize()


    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    eps = 1e-3
    # We then retrieve all the volumes in the bounding box of the original cube,
    # and delete all the parts outside it:
    vin = gmsh.model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 3)
    for v in vin
        deleteat!(out, findall(x -> x == v, out))
    end
    # Delete outside parts recursively
    gmsh.model.removeEntities(out, true)


    return gmsh.model
end

function getBoundingBox(rve::RVE)
    xmin = rve.origin[1]
    ymin = rve.origin[2]
    zmin = rve.origin[3]
    xmax = rve.origin[1] + rve.size[1]
    ymax = rve.origin[2] + rve.size[2]
    zmax = rve.origin[3] + rve.size[3]
    return [xmin, ymin, zmin, xmax, ymax, zmax]
end

function createMesh!(model::Module, rve::RVE, inclusions::Tuple,  out_::Vector{Any} )
    box = getBoundingBox(rve)
    eps = 1e-3
    # vin = model.getEntitiesInBoundingBox(box[1] - eps, box[2] - eps, box[3] - eps,
    #     box[4] + eps, box[5] + eps, box[6] + eps, 3)
    # p = model.getBoundary(vin, false, false, true)  # Get all points
    # model.mesh.setSize(p, meshsize)
    # p = model.getEntitiesInBoundingBox(box[1] - eps, box[2] - eps, box[3] - eps, box[1] + eps, box[2] + eps, box[3] + eps, 0)
    # model.mesh.setSize(p, 0.001)

 # mesh matrix
incu = 1
taginc = model.getEntitiesForPhysicalGroup(3, incu)
model.mesh.field.add("MathEval", 1)
model.mesh.field.setString(1, "F",
    "$(rve.meshsize)")
model.mesh.field.add("Restrict", 2)
model.mesh.field.setNumber(2, "InField", 1)
model.mesh.field.setNumbers(2, "VolumesList", vcat(taginc...))

numfield = 2
for i in 1:length(out_)-1
        taginc = model.getEntitiesForPhysicalGroup(3,  i+2)
        numfield += 1
        model.mesh.field.add("MathEval", numfield)
        model.mesh.field.setString(numfield, "F",
        "$(inclusions[i].meshsize)")
        numfield += 1
        model.mesh.field.add("Restrict", numfield)
        model.mesh.field.setNumber(numfield, "InField", numfield - 1)
        model.mesh.field.setNumbers(numfield, "VolumesList", vcat(taginc...))
end

numfield += 1
model.mesh.field.add("Min", numfield)
model.mesh.field.setNumbers(numfield, "FieldsList", collect(2:2:numfield-1))
model.mesh.field.setAsBackgroundMesh(numfield)


gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)




    sxmin = model.getEntitiesInBoundingBox(box[1] - eps, box[2] - eps, box[3] - eps,
        box[4] + eps, box[5] + eps, box[6] + eps, 2)

    if rve.periodicityFlags[1] == 1
        println("periodicity x!!")
        for i in sxmin
            # Then we get the bounding box of each left surface
            xmin, ymin, zmin, xmax, ymax, zmax = model.getBoundingBox(i[1], i[2])
            # We translate the bounding box to the right and look for surfaces inside
            # it:
            sxmax = model.getEntitiesInBoundingBox(xmin - eps + rve.size[1], ymin - eps,
                zmin - eps, xmax + eps + rve.size[1],
                ymax + eps, zmax + eps, 2)
            # For all the matches, we compare the corresponding bounding boxes...
            for j in sxmax
                xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = model.getBoundingBox(
                    j[1], j[2])
                xmin2 -= rve.size[1]
                xmax2 -= rve.size[1]
                # ...and if they match, we apply the periodicity constraint
                if (abs(xmin2 - xmin) < eps && abs(xmax2 - xmax) < eps &&
                    abs(ymin2 - ymin) < eps && abs(ymax2 - ymax) < eps &&
                    abs(zmin2 - zmin) < eps && abs(zmax2 - zmax) < eps)
                    model.mesh.setPeriodic(2, [j[2]], [i[2]], [1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
                end
            end
        end
    end

    if rve.periodicityFlags[2] == 1
        println("periodicity y!!")
        for i in sxmin
            # Then we get the bounding box of each left surface
            xmin, ymin, zmin, xmax, ymax, zmax = model.getBoundingBox(i[1], i[2])
            # We translate the bounding box to the right and look for surfaces inside
            # it:
            sxmax = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps + rve.size[2],
                zmin - eps, xmax + eps,
                ymax + eps + rve.size[2], zmax + eps, 2)
            # For all the matches, we compare the corresponding bounding boxes...
            for j in sxmax
                xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = model.getBoundingBox(
                    j[1], j[2])
                ymin2 -= rve.size[2]
                ymax2 -= rve.size[2]
                # ...and if they match, we apply the periodicity constraint
                if (abs(xmin2 - xmin) < eps && abs(xmax2 - xmax) < eps &&
                    abs(ymin2 - ymin) < eps && abs(ymax2 - ymax) < eps &&
                    abs(zmin2 - zmin) < eps && abs(zmax2 - zmax) < eps)
                    model.mesh.setPeriodic(2, [j[2]], [i[2]], [1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1])
                end
            end
        end
    end

    if rve.periodicityFlags[3] == 1
        println("periodicity z!!")
        for i in sxmin
            # Then we get the bounding box of each left surface
            xmin, ymin, zmin, xmax, ymax, zmax = model.getBoundingBox(i[1], i[2])
            # We translate the bounding box to the right and look for surfaces inside
            # it:
            sxmax = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps,
                zmin - eps + rve.size[3], xmax + eps,
                ymax + eps, zmax + eps + rve.size[3], 2)
            # For all the matches, we compare the corresponding bounding boxes...
            for j in sxmax
                xmin2, ymin2, zmin2, xmax2, ymax2, zmax2 = model.getBoundingBox(
                    j[1], j[2])
                zmin2 -= rve.size[3]
                zmax2 -= rve.size[3]
                # ...and if they match, we apply the periodicity constraint
                if (abs(xmin2 - xmin) < eps && abs(xmax2 - xmax) < eps &&
                    abs(ymin2 - ymin) < eps && abs(ymax2 - ymax) < eps &&
                    abs(zmin2 - zmin) < eps && abs(zmax2 - zmax) < eps)
                    model.mesh.setPeriodic(2, [j[2]], [i[2]], [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1])
                end
            end
        end
    end

    model.mesh.generate(3)

    return model
end

function visualizeMesh()
    # Launch the GUI to see the results:
    if !("-nopopup" in ARGS)
        gmsh.fltk.run()
    end

end

function closeGmshModel()
    Gmsh.finalize()
end


function saveMesh(model::Module, output_file::String)
    gmsh.write(output_file)
end



function _addPhysicalGroups!(model::Module, out_::Vector{Any}, rve::RVE)
    xmin, ymin, zmin, xmax, ymax, zmax = getBoundingBox(rve)
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




function _removeEntities!(model::Module, rve::RVE, out::Array{Tuple{Int32,Int32},1}, recursive::Bool)
    xmin, ymin, zmin, xmax, ymax, zmax = getBoundingBox(rve)
    gmsh.option.setNumber("Geometry.OCCBoundsUseStl", 1)
    eps = 1e-3
    vin = model.getEntitiesInBoundingBox(xmin - eps, ymin - eps, zmin - eps, xmax + eps, ymax + eps, zmax + eps, 3)
    for v in vin
        deleteat!(out, findall(x -> x == v, out))
    end
    # Delete outside parts recursively
    model.removeEntities(out, recursive)
    model.occ.synchronize()
    visualizeMesh()

end






rve = RVE([2, 2, 2], [1.0, 1.0, 1.0], [1, 1, 1], [0.0, 0.0, 0.0], 0.1)
inclusions = (Sphere([1.0, 1.0, 0.5], 0.25, 0.01), Sphere([0.5, 0.5, 1.0], 0.25, 0.05), Sphere([0.5, 0.5, 0.5], 0.2, 0.03))
# inclusions = (Sphere([1.0, 1.0, 0.5], 0.25, 0.1), Cylinder([0.5, 0.5, 0.1], [0.0, 0.0, 0.9], 0.25, 0.1))




gmsh.initialize()
gmsh.model.add("test")
xmin, ymin, zmin, xmax, ymax, zmax = getBoundingBox(rve)
_addBoundingBox!(gmsh.model, rve)
dim3::Int32 = 3
numinc = length(inclusions)
_tags = Vector{Vector{Tuple{Int32,Int32}}}()
push!(_tags, [(3, 1)])

_numinc = zeros(Int32, numinc)
for i in 1:numinc
    tag = _addInclusion!(gmsh.model, inclusions[i])
    # gmsh.model.occ.synchronize()
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
push!(out_, out__[1][2:end])
for i in 2:numinc
    _, out__ = gmsh.model.occ.fragment(out__[1][1], vcat(_tags[i+1]...))
    push!(out_, out__[1][2:end])
end
push!(out_, out__[1][1])
out = vcat(out_...)

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
  
createMesh!(gmsh.model, rve, inclusions,  out_ )

visualizeMesh()



 






 

