struct MeshStrategy{Kind} end

function createMesh!(model::Module, rve::RVE, inclusion::Inclusion, out_::Vector{Any})
    return createMesh!(model, rve, (inclusion,), out_)
end


function createMesh!(model::Module, rve::RVE, inclusions::Tuple, out_::Vector{Any})
    box = _getBoundingBox(rve)
    eps = 1e-3
 
    field_list=[]
    incu = 1
    taginc = model.getEntitiesForPhysicalGroup(3, incu)
    model.mesh.field.add("MathEval", 1)
    model.mesh.field.setString(1, "F",
        "$(rve.meshsize)")
    model.mesh.field.add("Restrict", 2)
    model.mesh.field.setNumber(2, "InField", 1)
    model.mesh.field.setNumbers(2, "VolumesList", vcat(taginc...))
    push!(field_list,2)
    
    numfield = 2
    for i in 1:length(out_)-1
        taginc = model.getEntitiesForPhysicalGroup(3, i + 2)
        tagdown = []
        for j in 1:length(taginc)
            _, tagdown_ = model.getAdjacencies(3, taginc[j])
            push!(tagdown, tagdown_)
        end

        lc=_getlc(inclusions[i])
        ref=_getreflevel(inclusions[i])  
        meshsize=2*pi*lc/15

        numfield += 1
        model.mesh.field.add("MathEval", numfield)
        model.mesh.field.setString(numfield, "F",
            "$(meshsize/ref)")
        numfield += 1
        model.mesh.field.add("Restrict", numfield)
        model.mesh.field.setNumber(numfield, "InField", numfield - 1)
        # model.mesh.field.setNumbers(numfield, "VolumesList", vcat(taginc...))
         model.mesh.field.setNumbers(numfield, "SurfacesList", vcat(tagdown...))
        push!(field_list,numfield)
    
        numfield += 1
        model.mesh.field.add("Distance", numfield)
        model.mesh.field.setNumbers(numfield, "SurfacesList", vcat(tagdown...))
        model.mesh.field.setNumber(numfield, "Sampling", 100)
    
        numfield += 1
        model.mesh.field.add("Threshold", numfield)
        model.mesh.field.setNumber(numfield, "InField", numfield-1)
        model.mesh.field.setNumber(numfield, "SizeMin", meshsize/ref)
        model.mesh.field.setNumber(numfield, "SizeMax", rve.meshsize)
        model.mesh.field.setNumber(numfield, "DistMin", lc / 5)
        model.mesh.field.setNumber(numfield, "DistMax", lc / 2.5)
        push!(field_list,numfield)
    end
    
    
    numfield += 1
    model.mesh.field.add("Min", numfield)
    model.mesh.field.setNumbers(numfield, "FieldsList", field_list)
    model.mesh.field.setAsBackgroundMesh(numfield)
    
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    


    sxmin = model.getEntitiesInBoundingBox(box[1] - eps, box[2] - eps, box[3] - eps,
        box[4] + eps, box[5] + eps, box[6] + eps, 2)

    if rve.periodicityFlags[1] == 1
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
 