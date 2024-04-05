"""
Set of Inclusions
"""
abstract type Inclusion end

struct Box <: Inclusion
    size::Vector{Float64}    # size of the RVE
    origin::Vector{Float64}     # set cell origin to [0,0,0]
    refinementwidth::Int32
end

struct Sphere <: Inclusion
    origin::Vector{Float64} # center of the sphere
    radius::Float64 # radius of the sphere
    refinementwidth::Int32
end

struct Cylinder <: Inclusion
    origin::Vector{Float64} # center of the sphere
    axis::Vector{Float64}  # axis of cilinder
    radius::Float64 #radius of cylinder
    refinementwidth::Int32
end

struct Ellipsoid <: Inclusion
    origin::Vector{Float64} # center of the sphere
    radius::Vector{Float64} # radius of the sphere
    θ::Vector{Float64} # [θxz, θxy]
    refinementwidth::Int32
end


struct Fuse <: Inclusion
    inclusions::Vector{Inclusion}
end

struct Cut <: Inclusion
    object::Inclusion
    tool::Inclusion
end

struct Intersect <: Inclusion
    object::Inclusion
    tool::Inclusion
end

struct Translate <: Inclusion
    inclusions::Inclusion
    dx::Float64
    dy::Float64
    dz::Float64
end


struct Rotate <: Inclusion
    inclusions::Inclusion
    origin::Vector{Float64}
    direction::Vector{Float64}
    ϕ::Float64
end

"""
Add inclusion to gmsh model
"""
function (inc::Translate)(model::Module)
    tag1_ = inc.inclusions(model)
    model.occ.translate((3, tag1_), inc.dx, inc.dy, inc.dz)
    return tag1_
end

function (inc::Rotate)(model::Module)
    tag1_ = inc.inclusions(model)
    model.occ.rotate((3, tag1_), inc.origin[1], inc.origin[2], inc.origin[3], inc.direction[1], inc.direction[2], inc.direction[3], inc.ϕ)
    return tag1_
end

function (inc::Fuse)(model::Module)
    tag1_ = inc.inclusions[1](model)
    for i in inc.inclusions[2:end]
        tag2_ = i(model)
        model.occ.fuse((3, tag1_), (3, tag2_), -1)
    end
    return tag1_
end

function (inc::Cut)(model::Module)
    tag1_ = inc.object(model)
    tag2_ = inc.tool(model)
    model.occ.cut((3, tag1_), (3, tag2_), -1)
    return tag1_
end

function (inc::Intersect)(model::Module)
    tag1_ = inc.object(model)
    tag2_ = inc.tool(model)
    model.occ.intersect((3, tag1_), (3, tag2_), -1)
    return tag1_
end

function (inc::Box)(model::Module)
    tag = model.occ.addBox(inc.origin[1], inc.origin[2], inc.origin[3], inc.size[1], inc.size[2], inc.size[3], -1)
    return tag
end

function (inc::Sphere)(model::Module)
    tag = model.occ.addSphere(inc.origin[1], inc.origin[2], inc.origin[3], inc.radius, -1)
    return tag
end

function (inc::Cylinder)(model::Module)
    tag = model.occ.addCylinder(inc.origin[1], inc.origin[2], inc.origin[3], inc.axis[1], inc.axis[2], inc.axis[3], inc.radius, -1)
    return tag
end

function (inc::Ellipsoid)(model::Module)
    tag = model.occ.addSphere(inc.origin[1], inc.origin[2], inc.origin[3], 1, -1)
    model.occ.dilate((3, tag), inc.origin[1], inc.origin[2], inc.origin[3], inc.radius[1], inc.radius[2], inc.radius[3])
    model.occ.rotate((3, tag), inc.origin[1], inc.origin[2], inc.origin[3], 0, 1, 0, inc.θ[1])
    model.occ.rotate((3, tag), inc.origin[1], inc.origin[2], inc.origin[3], 0, 0, 1, inc.θ[2])
    return tag
end



"""
Get characteristic length for meshing
"""

function _getlc(inc::Sphere)
    return inc.radius
end

function _getlc(inc::Cylinder)
    return inc.radius
end

function _getlc(inc::Ellipsoid)
    return minimum(inc.radius)
end

function _getlc(inc::Box)
    return minimum(inc.size)
end

function _getlc(inc::Fuse)
    radius = zeros(Float64, length(inc.inclusions))
    for i in 1:length(inc.inclusions)
        radius[i] = _getlc(inc.inclusions[i])
    end
    return minimum(radius)
end

function _getlc(inc::Cut)
    radius1 = _getlc(inc.object)
    radius2 = _getlc(inc.tool)
    return minimum([radius1, radius2])
end

function _getlc(inc::Intersect)
    radius1 = _getlc(inc.object)
    radius2 = _getlc(inc.tool)
    return minimum([radius1, radius2])
end

function _getlc(inc_::Translate)
    return _getlc(inc_.inclusions)
end

function _getlc(inc_::Rotate)
    return _getlc(inc_.inclusions)
end

"""
Get refinement level
"""


function _getreflevel(inc::Sphere)
    return inc.refinementwidth
end

function _getreflevel(inc::Cylinder)
    return inc.refinementwidth
end

function _getreflevel(inc::Ellipsoid)
    return inc.refinementwidth
end

function _getreflevel(inc::Box)
    return inc.refinementwidth
end

function _getreflevel(inc::Fuse)
    ref = zeros(Float64, length(inc.inclusions))
    for i in 1:length(inc.inclusions)
        ref[i] = _getreflevel(inc.inclusions[i])
    end
    return maximum(ref)
end

function _getreflevel(inc::Cut)
    ref1 = _getreflevel(inc.object)
    ref2 = _getreflevel(inc.tool)
    return minimum([ref1, ref2])
end

function _getreflevel(inc::Intersect)
    ref1 = _getreflevel(inc.object)
    ref2 = _getreflevel(inc.tool)
    return minimum([ref1, ref2])
end

function _getreflevel(inc_::Translate)
    return _getreflevel(inc_.inclusions)
end

function _getreflevel(inc_::Rotate)
    return _getreflevel(inc_.inclusions)
end