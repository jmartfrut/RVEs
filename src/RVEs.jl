module RVEs

using Gmsh: Gmsh, gmsh

export RVE
export Sphere
export Cylinder
export Ellipsoid
export Box
export Fuse
export Cut
export Intersect 
export Translate
export Rotate
export createGmshModel
export stopGmsh
export startGmsh
export visualizeMesh
export createMesh!
export saveMesh
export ShowInfo
export MeshStrategy

struct RVE
    size::Vector{Float64}    # size of the RVE
    periodicityFlags::Vector{Int64}       # periodicity flags
    origin::Vector{Float64}     # set cell origin to [0,0,0]
    meshsize::Float64
end
 

include("Gmsh_Inclusions.jl")
include("Gmsh_Utils.jl")
include("Gmsh_Mesh.jl")
include("Gmsh_Model.jl")
 
 

end
