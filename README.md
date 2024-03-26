# RVEs :construction: :construction: :construction: **Work in progress** :construction: :construction: :construction:


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jmartfrut.github.io/RVEs.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jmartfrut.github.io/RVEs.jl/dev/)
[![Build Status](https://github.com/jmartfrut/RVEs.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jmartfrut/RVEs.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/jmartfrut/RVEs.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jmartfrut/RVEs.jl)

# Modelling of periodic RVEs using GMSH

**Gmsh:**
Gmsh stands as an open-source 3D finite element mesh generator equipped with an integrated CAD engine and post-processor. It's designed for speed, efficiency, and user-friendliness, featuring parametric input and adaptable visualization capabilities. The software architecture revolves around four core modules: geometry, mesh, solver, and post-processing. These modules offer flexible control options, including a graphical interface, command-line interface, dedicated .geo files, and support for APIs in C++, C, Python, Julia, and Fortran.

**RVEs Package:**
The RVEs package serves as a versatile tool for mesh generation, leveraging the Gmsh-Julia-API for seamless integration. It provides a robust framework for automating mesh generation for custom model types, enabling users to automate the creation of intricate structures such as representative volume elements (RVEs). RVEs simplifies the mesh modeling process into several key stages:

1. **Geometry Setup:** Define geometric shapes using fundamental entities and boolean operations.
2. **Integration with Gmsh:** Incorporate shapes into Gmsh, execute boolean operations, and define physical groups.
3. **Mesh Generation:** Generate meshes with customizable refinement fields.
4. **Export and Visualization:** Export meshes into various formats and visualize the results.

**Current Focus and Future Expansion:**
At present, RVEs excels in automating the generation of representative volume elements containing multiple inclusion objects. However, its capabilities can be extended by adding support for additional geometric shapes and model types. It's essential to clarify that RVEs does not aim to replace the Gmsh scripting language. Instead, it acts as a facilitator, simplifying the automation of complex meshing tasks within a user-friendly programming environment such as Julia.

# MiniApps included  

 <p align="center"> 
&nbsp; &nbsp; &nbsp; &nbsp;
<img alt="Dark"
src="https://github.com/jmartfrut/RVEs/blob/main/docs/imgs/rves.png" width="100%">
</p>






# Quick start

```julia
using RVEs

startGmsh()

myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)

myInclusions = (Sphere([0.0, 0.5, 0.5], 0.25, 2),
Sphere([0.5, 0.0, 0.5], 0.25, 2),
Sphere([0.5, 0.5, 0.0], 0.25, 2),
Sphere([1.0, 1.0, 1.0], 0.25, 2))

model, tagvol=createGmshModel(myRVE,myInclusions,"RVE")

createMesh!(model,myRVE,myInclusions,tagvol)

output_file = joinpath(dirname(@__FILE__), "RVE.msh")
saveMesh(output_file)

stopGmsh()
```



```julia
using RVEs

startGmsh()

myRVE = RVE( [1.0, 1.0, 1.0], [1,1,1], [0.0, 0.0, 0.0], 0.1)
 
geo1=Cylinder([-0.1, 0.5, 0.5],[2.0,0.0,0.0], 0.1, 2)
geo2=Cylinder([0.5, -0.1, 0.5],[0.0,2.0,0.0], 0.1, 2)
geo3=Cylinder([0.5, 0.5,-0.1],[0.0,0.0,2.0], 0.1, 2)
geo4=Fuse([geo1,geo2,geo3],2)

geo5=Sphere([0.5, 0.5, 0.5], 0.25, 2)
L=0.38
geo6=Box([L, L, L], [0.5-L/2, 0.5-L/2, 0.5-L/2], 2)
geo7=Intersect(geo5,geo6,2)

geo8=Cut(geo7,geo4,2)
myInclusions=(geo8,)
  
model, tagvol=createGmshModel(myRVE,myInclusions,"RVE")
  
createMesh!(model,myRVE,myInclusions,tagvol)

output_file = joinpath(dirname(@__FILE__), "RVE.msh")
saveMesh(output_file)

stopGmsh()
```


# Project funded by:

- Grant PID2022-141957OA-C22 funded by MCIU/AEI/ 10.13039/501100011033  and by ''ERDF A way of making Europe''

 <p align="center"> 
&nbsp; &nbsp; &nbsp; &nbsp;
<img alt="Dark"
src="https://github.com/jmartfrut/RVEs/blob/main/docs/imgs/aei.png" width="60%">
</p>