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
