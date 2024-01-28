module Rect_waveguide_quads

# Guia de onda rectangular, metodo de elementos finitos con cuadrilateros

import gmsh

gmsh.initialize()

gmsh.model.add("guiaonda_rect_cuadrilateros")

tm = 1.0
N = 8

# Definir puntos
p1 = gmsh.model.geo.addPoint(0, 0, 0, tm)
p2 = gmsh.model.geo.addPoint(10, 0, 0, tm)
p3 = gmsh.model.geo.addPoint(10, 5, 0, tm)
p4 = gmsh.model.geo.addPoint(0, 5, 0, tm)

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
gmsh.model.addPhysicalGroup(1, [1, 2, 3, 4], 5, "DIRICHLET")
#gmsh.model.addPhysicalGroup(1, [1, 2], 5, "DIRICHLET")
#gmsh.model.addPhysicalGroup(1, [3, 4], 6, "NEUMAN")
#gmsh.model.addPhysicalGroup(1, [1, 2], 5, "DIRICHLET")
#gmsh.model.addPhysicalGroup(1, [4], 6, "NEUMAN")
#gmsh.model.addPhysicalGroup(1, [3], 7, "CAUCHY")
gmsh.model.addPhysicalGroup(2, [1], 6, "Aire dentro de la guia")

# Configura el mallado y la conversion a cuadrilateros
#gmsh.model.mesh.setTransfiniteCurve(1, tm+1)  # El número de nodos será N+1
#gmsh.model.mesh.setTransfiniteSurface(1)
gmsh.model.mesh.generate(2)
gmsh.model.mesh.recombine()

# Generar los elementos de orden 2
#gmsh.model.mesh.setOrder(2)

# Obtenemos los identificadores de los nodos, las coordenadas y las coordenadas parametricas de los mismos
export nodeTags_c, coord_c
nodeTags_c, coord_c = gmsh.model.mesh.getNodes()

# Obtenemos los tipos de elementos, sus identificadores y los nodos de todos los elementos del mismo tipo concatenados,
# estos son los identificadores unicos de los nodos que forman parte de los elementos.
export elementTypes_c, elementTags_c, nodeTags_e_c, element_prop_c
elementTypes_c, elementTags_c, nodeTags_e_c = gmsh.model.mesh.getElements()
element_prop_c = gmsh.model.mesh.getElementProperties(elementTypes_c[2])

# Obtenemos las entidades
export entities_c, physical_group_c, boundary_names_c, physical_names_c, physical_material_c, materials_names_c
entities_c = gmsh.model.getEntities()

physical_group_c = gmsh.model.getPhysicalGroups(1)
boundary_names_c = Vector{Tuple{Int32, String}}()
for group in physical_group_c
    boundary_name = gmsh.model.getPhysicalName(group[1], group[2])
    push!(boundary_names_c, (group[2], boundary_name))
end

physical_material_c = gmsh.model.getPhysicalGroups(2)
materials_names_c = Vector{Tuple{Int32, String}}()
for material in physical_material_c 
    material_name = gmsh.model.getPhysicalName(material[1], material[2])
    push!(materials_names_c, (material[2], material_name))
end

# Obtenemos los nodos que pertenecen a los grupos fisicos
export groups_c
groups_c = Vector{Tuple{Int32, Vector{UInt64}}}()
for group in physical_group_c
    phy_nodes = gmsh.model.mesh.getNodesForPhysicalGroup(1, group[2])
    push!(groups_c, (group[2], phy_nodes[1]))
end

# Finalizar el modelo gmsh
gmsh.finalize()

end