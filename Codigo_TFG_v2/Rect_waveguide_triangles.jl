module Rect_waveguide_triangles

# Guia de onda rectangular metodo de elementos finitos con triangulos

import gmsh

gmsh.initialize()

gmsh.model.add("rect_waveguide_triangles")

#=tm = 1.0
global a = 10
global b = 5=#

global a = 0.02286
global b = 0.01016
tm = a/5

#= Definir puntos
p1 = gmsh.model.geo.addPoint(0, 0, 0, tm)
p2 = gmsh.model.geo.addPoint(a, 0, 0, tm)
p3 = gmsh.model.geo.addPoint(a, b, 0, tm)
p4 = gmsh.model.geo.addPoint(0, b, 0, tm)=#
p1 = gmsh.model.geo.addPoint(-a/2, -b/2, 0, tm)
p2 = gmsh.model.geo.addPoint(a/2, -b/2, 0, tm)
p3 = gmsh.model.geo.addPoint(a/2, b/2, 0, tm)
p4 = gmsh.model.geo.addPoint(-a/2, b/2, 0, tm)

# Definir lineas
gmsh.model.geo.addLine(p1, p2, 1)
gmsh.model.geo.addLine(p2, p3, 2)
gmsh.model.geo.addLine(p3, p4, 3)
gmsh.model.geo.addLine(p4, p1, 4)

# Definir bucles de curvas
gmsh.model.geo.addCurveLoop([1, 2, 3, 4], 1)

# Definir superficies
gmsh.model.geo.addPlaneSurface([1], 1)

gmsh.model.geo.synchronize()

# Definir grupos fisicos

# Solo una condicion de contorno
gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4], 5, "DIRICHLET")

# Dos condiciones de contorno 
#gmsh.model.addPhysicalGroup(1, [1, 2], 5, "DIRICHLET")
#gmsh.model.addPhysicalGroup(1, [3, 4], 6, "NEUMAN")

# Tres condiciones de contorno
#gmsh.model.addPhysicalGroup(1, [1, 2], 5, "DIRICHLET")
#gmsh.model.addPhysicalGroup(1, [4], 6, "NEUMAN")
#gmsh.model.addPhysicalGroup(1, [3], 7, "CAUCHY")

gmsh.model.addPhysicalGroup(2, [1], 6, "Aire dentro de la guia")

# Generar la malla
gmsh.model.mesh.generate(2)

# Generar los elementos de orden 2
#gmsh.model.mesh.setOrder(2)

# Obtenemos los identificadores de los nodos, las coordenadas
export nodeTags_t, coord_t
nodeTags_t, coord_t = gmsh.model.mesh.getNodes()

#= Obtenemos los tipos de elementos, sus identificadores 
y los nodos de todos los elementos del mismo tipo concatenados,
estos son los identificadores unicos de los nodos que forman parte de los elementos.=#
export elementTypes_t, elementTags_t, nodeTags_e_t, element_prop_t
elementTypes_t, elementTags_t, nodeTags_e_t = gmsh.model.mesh.getElements()
element_prop_t = gmsh.model.mesh.getElementProperties(elementTypes_t[2])

# Obtenemos las entidades
export entities_t, physical_group_t, boundary_names_t, physical_names_t, physical_material_t, materials_names_t
entities_t = gmsh.model.getEntities()
physical_group_t = gmsh.model.getPhysicalGroups(1)

boundary_names_t = Vector{Tuple{Int32, String}}()
for group in physical_group_t
    boundary_name = gmsh.model.getPhysicalName(group[1], group[2])
    push!(boundary_names_t, (group[2], boundary_name))
end

physical_material_t = gmsh.model.getPhysicalGroups(2)
materials_names_t = Vector{Tuple{Int32, String}}()
for material in physical_material_t 
    material_name = gmsh.model.getPhysicalName(material[1], material[2])
    push!(materials_names_t, (material[2], material_name))
end

# Obtenemos los nodos que pertenecen a los grupos fisicos
export groups_t
groups_t = Vector{Tuple{Int32, Vector{UInt64}}}()
for group in physical_group_t
    phy_nodes = gmsh.model.mesh.getNodesForPhysicalGroup(1, group[2])
    push!(groups_t, (group[2], phy_nodes[1]))
end

# Finalizar el modelo gmsh
gmsh.finalize()

end