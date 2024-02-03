module Mesh_object

include("Rect_waveguide_triangles.jl")
include("Rect_waveguide_quads.jl")
include("Elements.jl")
include("Elements_order.jl")
include("Material.jl")
include("Boundary_cond.jl")

# Variable para cambiar de mallas triangulares a mallas cuadradas
###############CAMBIAR####################
quadrangular_meshes = false;

# Vector con todas las condiciones de contorno de la malla
mesh_cond = Vector{Boundary_cond.Boundary_condition}()

# Vector con todos los materiales de la malla
mesh_materials = Vector{Material.My_material}()

# Vector con todas las coordenadas de los nodos de la malla
mesh_coord_2D = Vector{Float64}()

# Vector con todos los elementos de la malla
mesh_elements = Vector{Elements.Element}()

# Funcion que obtiene todos los materiales de la malla
function get_all_materials(materials_names::Vector{Tuple{Int32, String}})
    """
        get_all_materials(materials_names)

    Obtiene todos los materiales de la malla 
    """
    for materials in materials_names
        material = Material.My_material(materials[1], materials[2], 0.0, 0.0)
        push!(mesh_materials, material)
    end
end

###############CAMBIAR####################
if quadrangular_meshes == false
    get_all_materials(Rect_waveguide_triangles.materials_names_t) # Clase geometria de entrada
else
    get_all_materials(Rect_waveguide_quads.materials_names_c)
end

# Funcion que obtiene todas las condiciones de contorno
function get_all_boundary_cond(boundary_names::Vector{Tuple{Int32, String}}, groups::Vector{Tuple{Int32, Vector{UInt64}}})
     """
        get_all_boundary_cond(boundary_names, groups)

    Obtiene todas las condiciones de contorno de la malla, asi como los nodos que las componen 
    """
    for boundarys in boundary_names
        nodes = Vector{UInt64}()
        for group in groups
            if group[1] == boundarys[1]
                nodes = group[2]
            end
        end
        if boundarys[2] == "DIRICHLET"
            type = Boundary_cond.DIRICHLET
        elseif boundarys[2] == "NEUMAN"
            type = Boundary_cond.NEUMAN
        end
    boundary = Boundary_cond.Boundary_condition(boundarys[1], type, nodes)
    push!(mesh_cond, boundary)
    end
end

###############CAMBIAR####################
if quadrangular_meshes == false
    get_all_boundary_cond(Rect_waveguide_triangles.boundary_names_t, Rect_waveguide_triangles.groups_t)
else
    get_all_boundary_cond(Rect_waveguide_quads.boundary_names_c, Rect_waveguide_quads.groups_c)
end

# Funcion que guarda las coordenadas x e y de todos los nodos de la malla en mesh_coord_2D
function get_all_coordenates(coord::Vector{Float64})
     """
        get_all_coordenates(coord)

    Guarda las coordenadas 2D de todos los nodos de la malla con la siguiente estructura: {n1x, n1y, n2x, n2y, ...} 
    """

    # Coord contiene las coordenadas completas con {n1x, n1y, n1z, n2x, n2y, n2z, ...}

    # Iteramos a traves de las coordenadas
    for i in 1:3:length(coord)  # Saltamos de 3 en 3 para seleccionar x, y
        x = coord[i]
        y = coord[i + 1]

        push!(mesh_coord_2D, x)
        push!(mesh_coord_2D, y)
    end
end

###############CAMBIAR####################
if quadrangular_meshes == false
    get_all_coordenates(Rect_waveguide_triangles.coord_t)
else
    get_all_coordenates(Rect_waveguide_quads.coord_c)
end

# Funcion para determinar el numero de nodos segun el orden del elemento
function evaluate_order(orden::Int32)
    """
        evaluate_order(orden)

    Determina el número de nodos que tiene el elemento según su orden 
    """        
    switch_dict = Dict(
        2 => 3,
        9 => 6,
        21 => 10,
        3 => 4,
        10 => 9
    )

    if haskey(switch_dict, orden)
        return switch_dict[orden]
    else
        throw(ArgumentError("No existe el orden definido"))
    end
end

# Funcion para determinar que nodos pertenecen a las condiciones de contorno
function get_boundary_nodes(element::Elements.Element, boundary_conditions::Vector{Boundary_cond.Boundary_condition})
    """
        get_boundary_nodes(element, boundary_conditions)

    Determina que nodos dentro del argumento "element" cumplen las condiciones de contorno  
    """
    if isa(element, Elements.My_triangle)
        edge_1 = [element.nodes[1], element.nodes[2]]
        edge_2 = [element.nodes[1], element.nodes[3]]
        edge_3 = [element.nodes[2], element.nodes[3]]
        edges = (edge_1, edge_2, edge_3)
    elseif isa(element, Elements.My_quad)  
        edge_1 = [element.nodes[1], element.nodes[2]]
        edge_2 = [element.nodes[2], element.nodes[3]]
        edge_3 = [element.nodes[3], element.nodes[4]]
        edge_4 = [element.nodes[1], element.nodes[4]]
        edges = (edge_1, edge_2, edge_3, edge_4)
    end
    
    boundary_nodes = Vector{Int64}()
    type_cond = 0
    for edge in edges
        for condition in boundary_conditions
            for i in edge
                for j in condition.nodes
                    if i==j
                        push!(boundary_nodes, i)
                    end
                end
            end
            if edge == boundary_nodes
                type_cond = condition.type
            end
            boundary_nodes =  Vector{Int64}()
        end
        if type_cond != 0
            push!(element.boundary_egdes, type_cond)
            type_cond = 0
        else
            push!(element.boundary_egdes, 0)
        end
    end
end

# Funcion que guarda todos los elementos de la malla en mesh_elements
function get_all_elements(elementType::Vector{Int32}, elemtTags::Vector{Vector{UInt64}}, nodeTags::Vector{Vector{UInt64}})
    """
        get_all_elements(elementType, elemtTags, nodeTags)

    Crea y guarda todos los elementos en el atributo "mesh_elements" junto con toda su información. 
    """
    for i in eachindex(elementType)
        if elementType[i] == Elements_order.GMSH_TRIA_ORDER_1 || elementType[i] == Elements_order.GMSH_TRIA_ORDER_2 || elementType[i] == Elements_order.GMSH_TRIA_ORDER_3
            # Numero de nodos por elemento 
            numNodesPerElement = evaluate_order(elementType[i])
            # Recorremos los elementos y creamos instancias de My_triangle
            elemtTags_e = elemtTags[i]
            for j in eachindex(elemtTags_e)
            # Calculamos el indice de inicio y final para los nodos de este elemento
            startIdx = (j - 1) * numNodesPerElement + 1
            endIdx = startIdx + numNodesPerElement - 1
        
            # Extraemos los nodos de este elemento
            nodeTags_e = nodeTags[i]
            nodesOfElement = nodeTags_e[startIdx:endIdx]
        
            # Creamos una instancia de My_triangle y la guardamos en elementObjects
            element = Elements.My_triangle(elemtTags_e[j],elementType[i], nodesOfElement, 0, Vector{Int64}(), Vector{Int64}())

            # Extraemos las condiciones de contorno
            get_boundary_nodes(element, mesh_cond)

            # Guardamos el elemento en mesh_elements
            push!(mesh_elements, element)
            end
        elseif elementType[i] == Elements_order.GMSH_QUAD_ORDER_1 || elementType[i] == Elements_order.GMSH_QUAD_ORDER_2
            # Numero de nodos por elemento 
            numNodesPerElement = evaluate_order(elementType[i])
            # Recorremos los elementos y creamos instancias de My_triangle
            elemtTags_e = elemtTags[i]
            for j in eachindex(elemtTags_e)
            # Calculamos el indice de inicio y final para los nodos de este elemento
            startIdx = (j - 1) * numNodesPerElement + 1
            endIdx = startIdx + numNodesPerElement - 1
        
            # Extraemos los nodos de este elemento
            nodeTags_e = nodeTags[i]
            nodesOfElement = nodeTags_e[startIdx:endIdx]

            # Creamos una instancia de My_triangle y la guardamos en elementObjects
            element = Elements.My_quad(elemtTags_e[j],elementType[i], nodesOfElement, 0, Vector{Int64}(), Vector{Int64}())

            # Extraemos las condiciones de contorno
            get_boundary_nodes(element, mesh_cond)

            # Guardamos el elemento en mesh_elements
            push!(mesh_elements, element)
            end
        end
    end
end

###############CAMBIAR####################
if quadrangular_meshes == false
    get_all_elements(Rect_waveguide_triangles.elementTypes_t, Rect_waveguide_triangles.elementTags_t, Rect_waveguide_triangles.nodeTags_e_t)
else
    get_all_elements(Rect_waveguide_quads.elementTypes_c, Rect_waveguide_quads.elementTags_c, Rect_waveguide_quads.nodeTags_e_c)
end

# Funcion que guarda los dof en cada elemento
function get_all_dofs(elements::Vector{Elements.Element}) 
    """
        get_all_dofs(elements)

    Asigna los grados de libertad a todos los nodos de la malla y devuelve el número de grados de libertad (dofs). 
    """
    shared_edges = Dict{Vector{Int64}, Int64}()
    shared_nodes = Dict{Int64, Int64}()
    cont = 0

    for element in elements

        nodes = element.nodes[1:3]
        for node in nodes
            if haskey(shared_nodes, node) 
                push!(element.dof, shared_nodes[node])
            else
                cont += 1
                push!(element.dof, cont)
                shared_nodes[node] = cont
            end
        end

        if (isa(element, Elements.My_triangle) && (element.order != 1))
            if isa(element, Elements.My_triangle)
                edge_1 = [element.nodes[1], element.nodes[2]]
                edge_2 = [element.nodes[1], element.nodes[3]]
                edge_3 = [element.nodes[2], element.nodes[3]]
                edges = (edge_1, edge_2, edge_3)

            elseif isa(element, Elements.My_quad)
                edge_1 = [element.nodes[1], element.nodes[2]]
                edge_2 = [element.nodes[2], element.nodes[3]]
                edge_3 = [element.nodes[3], element.nodes[4]]
                edge_4 = [element.nodes[1], element.nodes[4]]
                edges = (edge_1, edge_2, edge_3, edge_4)
            end
            for edge in edges
                if haskey(shared_edges, edge) 
                    push!(element.dof, shared_edges[edge])
                    if element.order == 3
                        push!(element.dof, shared_edges[edge]+1)
                    end
                else
                    cont += 1
                    push!(element.dof, cont)
                    shared_edges[edge] = cont
                    if element.order == 3
                        cont += 1
                        push!(element.dof, cont)
                        #shared_edges[edge] = cont
                    end
                end
            end
        end
         
        if (isa(element, Elements.My_triangle) && (element.order == 3)) || (isa(element, Elements.My_quad) && (element.order == 2))
            cont +=1
            push!(element.dof, cont)
        end
    end 
    return cont
end

total_dofs = get_all_dofs(mesh_elements)

end
