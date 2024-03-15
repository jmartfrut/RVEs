module MimosaRVE


using Gmsh: Gmsh, gmsh

export RVE
export Sphere
export Inclusion
export createGmshModel
export closeGmshModel
export visualizeMesh
export createMesh!
export saveMesh

abstract type Inclusion end

struct RVE
    nCells::Vector{Int64}    # number of cells in each direction
    size::Vector{Float64}    # size of the RVE
    periodicityFlags::Vector{Int64}       # periodicity flags
    origin::Vector{Float64}     # set cell origin to [0,0,0]
end

struct Sphere <: Inclusion
    origin::Vector{Float64} # center of the sphere
    radius::Float64 # radius of the sphere
end




function createGmshModel(rve::RVE, sf::Sphere, namedModel::String)
    gmsh.initialize()
    gmsh.model.add(namedModel)
    # coordinates of bounding box
    xmin = rve.origin[1]
    ymin = rve.origin[2]
    zmin = rve.origin[3]
    xmax = rve.origin[1] + rve.size[1]
    ymax = rve.origin[2] + rve.size[2]
    zmax = rve.origin[3] + rve.size[3]
    # create box. last number ID of entity
    gmsh.model.occ.addBox(xmin, ymin, zmin, xmax, ymax, zmax, 1)
    # create sphere. last number ID of entity
    gmsh.model.occ.addSphere(sf.origin[1], sf.origin[2], sf.origin[3], sf.radius, 2)

    #create periodic inclusions
    x̂ = -(xmin + xmax) * 0.5
    ŷ = -(ymin + ymax) * 0.5
    ẑ = -(zmin + zmax) * 0.5
    _tags = [(3, 2)]
    if rve.periodicityFlags[1] == 1
        newtags_x = gmsh.model.occ.copy(_tags)
        gmsh.model.occ.mirror(newtags_x, 1, 0, 0, x̂)
        _tags = vcat(_tags, newtags_x)
    end

    if rve.periodicityFlags[2] == 1
        newtags_y = gmsh.model.occ.copy(_tags)
        gmsh.model.occ.mirror(newtags_y, 0, 1, 0, ŷ)
        _tags = vcat(_tags, newtags_y)
    end

    if rve.periodicityFlags[2] == 1
        newtags_z = gmsh.model.occ.copy(_tags)
        gmsh.model.occ.mirror(newtags_z, 0, 0, 1, ẑ)
        _tags = vcat(_tags, newtags_z)
    end

    # fragment all the volumes, which will leave parts of spheres protruding outside the cube. second number IDs
    out, _ = gmsh.model.occ.fragment([(3, 1)], _tags)
    gmsh.model.occ.synchronize()
    # Ask OpenCASCADE to compute more accurate bounding boxes of entities using the STL mesh:    
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

function createMesh!(model::Module, rve::RVE, meshsize::Float64)
    box = getBoundingBox(rve)
    eps = 1e-3
    vin = model.getEntitiesInBoundingBox(box[1] - eps, box[2] - eps, box[3] - eps,
        box[4] + eps, box[5] + eps, box[6] + eps, 3)
    p = model.getBoundary(vin, false, false, true)  # Get all points
    model.mesh.setSize(p, meshsize)
    # p = model.getEntitiesInBoundingBox(box[1] - eps, box[2] - eps, box[3] - eps, box[1] + eps, box[2] + eps, box[3] + eps, 0)
    # model.mesh.setSize(p, 0.001)


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
            sxmax = model.getEntitiesInBoundingBox(xmin - eps , ymin - eps + rve.size[2],
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
            sxmax = model.getEntitiesInBoundingBox(xmin - eps , ymin - eps,
                zmin - eps+ rve.size[3], xmax + eps,
                ymax + eps , zmax + eps+ rve.size[3], 2)
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

# createMesh()

# saveMesh()

# visualizeMesh()









end
