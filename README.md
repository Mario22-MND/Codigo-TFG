# Code for 2D finite element simulator in Julia

This repository shows the code created to build a 2D finite element simulator in Julia. The most important information about the code will be documented below.

# Table of Contents

- [Installation](#installation)
- [Use](#use)
- [Future expansions](#futureExpansions)

# Installation
In order to use the code you must install several Julia packages and import the Gmsh API. The packages or libraries to install are: 
- LinearAlgebra: To obtain the eigenvalues and eigenvectors of the mass and stiffness matrix.
- SimplexQuad (https://github.com/eschnett/SimplexQuad.jl): allows us to obtain the integration points and the pessos according to the Gauss-Radau quadrature.
- Plots: to represent the meshes as well as the final results.
- LaTeXStrings: to mathematically customize our plots.

## Gmsh API
In order to use the Gmsh API we need the dynamic library Gmsh and the Julia module ('gmsh.jl'). This is obtained from the binary software development kit (SDK) available from the Gmsh website (https://gmsh.info/bin/), for Windows, Linux, and macOS. By downloading and unzipping the gmsh*-sdk.* file corresponding to your operating system and adding the directory (lib) directory of the SDK to JULIA_LOAD_PATH.  Once we have downloaded and unzipped the SDK, we look for the /.julia/config folder.
In case it does not exist, we will have to create it. Once the folder is created, create or edit the startup.jl file as follows:

```Julia
 # We save the path to the file we want to load in
JULIA_LOAD_PAD push !( LOAD_PATH , "/path/to/file ")
```

This startup.jl file is used to make custom configurations and load specific packages or modules that you want to be available every time you start Julia. The last step is to restart Julia.


# Use
The code consists of a series of modules, where the main one is called "Numerical_integration_v1.jl". This is responsible for calling "Mesh_object.jl" which in turn calls all the others. The dependencies of the modules are shown in the following image:

![dependencias_código](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/2b5e0d2e-07d6-4e69-8abd-a097c1722f26)

Where the arrows show from whom the information is obtained or on whom it depends.

## Rect_waveguide_triangles/quads.jl 
These two code snippets use the Gmsh API in Julia for the creation of the object to be used. From them we get important information such as nodes, node coordinates, physical groups (boundary conditions), etc.

## Elements.jl
This module is used to store the different types of elements. Both triangles and quads have been implemented. However, the orientation of quads has not been done, which is pending for future work. In this module "Elements_order.jl" is called where the order of the elements is stored transforming it to a readable format from the internal Gmsh codes shown in this web page: https://docs.juliahub.com/GmshTools/9rYp5/0.4.2/element_types/ .

## Material.jl
Thanks to this module, the different materials of which the mesh is composed can be saved. It will be useful for future work where generic meshes are analyzed. For the moment it has no influence on the code.

## Boundary_cond.jl
In this module we store the boundary conditions present in the mesh as well as all the nodes that suffer it. This module is used to check which nodes of the element have a certain boundary condition.

## Mesh_object.jl
Here you get all the important information to assemble the mesh. The elements are created where the boundary conditions are taken into account, the degrees of freedom are assigned, ect. The materials of which the mesh is composed are also saved as well as all the coordinates of the nodes that compose it.

All this information will be necessary to perform the calculations in "Numerical_integration_v1.jl".

### Mesh_plot.jl
In this module a graph of the assembled mesh is created using all the information created and obtained in "Mesh_object.jl". Examples of some meshes, of order 1 and order 2, are shown below:

![triangular_mesh_order_1_boundary_1](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/4d66a5eb-a339-4931-8c25-b615d9f5d46b)

![triangular_mesh_order_2_boundary_1](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/2c6c1603-00d9-4976-86e7-a3ff2871cfb7)

## Numerical_integration_v1.jl
This is the main module, since it is where all the calculations are performed with the information obtained so far and where the final results are obtained. It calculates the mass and stiffness matrices used to obtain the eigenvalues. You can choose by means of the attribute 'FEM_TE' if you want to calculate TE modes or TM modes. 

With the obtained eigenvalues, we calculate the relative error, which decreases as we increase the order of the elements and we make the mesh finer with the 'tm' attribute of the "Rect_waveguide_triangles.jl" module. We must take into account two points when calculating the relative error depending on whether we calculate TE or TM modes:
- TE modes: the first value obtained in 'k_c_fem' is used as reference, the second is the one corresponding to TE_10 mode.
  ![resultado_modos_TE](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/b73a9988-cc92-4e41-baee-d568cc5c5d65)

- TM modes: there are no TM_10 or TM_01 modes, so the first value obtained in 'k_c_fem' is the one corresponding to TM_11 mode.
  ![resultado_modos_TM](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/eaf11dc1-6092-47f8-9fae-81f0b001d839)

### Plot_results.jl
This module is used to plot the results and perform checks on them. Here we have checked that the error convergence rate is O(h^2p).

![Mode_TE](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/159bb6d6-a2c8-4664-ac10-9a783ff31b60)

![Mode_TM](https://github.com/Mario22-MND/Codigo-TFG/assets/126000794/b3a7919a-edba-4cfd-aba0-0469e59ad6ed)

# Future Expansions
The code can be extended by implementing the parts of the code dedicated to quadrilaterals, as well as implementing the use of materials to solve generic guides.
