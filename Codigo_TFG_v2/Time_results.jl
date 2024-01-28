

import gmsh

using LinearAlgebra
using SimplexQuad
using TimerOutputs

global GMSH_TRIA_ORDER_1 = 2
        global GMSH_TRIA_ORDER_2 = 9
        global GMSH_TRIA_ORDER_3 = 21
        global GMSH_QUAD_ORDER_1 = 3
        global GMSH_QUAD_ORDER_2 = 10

        abstract type Element end

        struct My_triangle <: Element
            id::Int64
            order::Int32
            nodes::Vector{Int64}
            material::Int64
            boundary_egdes::Vector{Int64}
            dof::Vector{Int64}

            # Constructor personalizado para Triangle
            function My_triangle(id::UInt64, order::Int32, nodes::Vector{UInt64}, material::Int64, boundary_egdes::Vector{Int64}, dof::Vector{Int64})
                # Determinar el orden de los elementos
                if order == GMSH_TRIA_ORDER_1
                    order = 1
                elseif order == GMSH_TRIA_ORDER_2 
                    order = 2
                elseif order == GMSH_TRIA_ORDER_3
                    order = 3
                else
                    throw(ArgumentError("El triangulo debe tener orden 1, 2 o 3"))
                end

                # Determinar las aristas de los elementos
                ord_nodes = sort(nodes[1:3])
                if order == 1
                    nodes = ord_nodes
                elseif order == 2
                    if ord_nodes[1] == nodes[1] && ord_nodes[2] == nodes[2] && ord_nodes[3] == nodes[3]
                        new_nodes = vcat(ord_nodes, nodes[4], nodes[6], nodes[5])
                    elseif ord_nodes[1] == nodes[1] && ord_nodes[2] == nodes[3] && ord_nodes[3] == nodes[2]
                        new_nodes = vcat(ord_nodes, nodes[6], nodes[4], nodes[5])
                    elseif ord_nodes[1] == nodes[2] && ord_nodes[2] == nodes[1] && ord_nodes[3] == nodes[3]
                        new_nodes = vcat(ord_nodes, nodes[4], nodes[5], nodes[6])
                    elseif ord_nodes[1] == nodes[2] && ord_nodes[2] == nodes[3] && ord_nodes[3] == nodes[1]
                        new_nodes = vcat(ord_nodes, nodes[5], nodes[4], nodes[6])
                    elseif ord_nodes[1] == nodes[3] && ord_nodes[2] == nodes[1] && ord_nodes[3] == nodes[2]
                        new_nodes = vcat(ord_nodes, nodes[6], nodes[5], nodes[4])
                    elseif ord_nodes[1] == nodes[3] && ord_nodes[2] == nodes[2] && ord_nodes[3] == nodes[1]
                        new_nodes = vcat(ord_nodes, nodes[5], nodes[6], nodes[4])
                    end
                    nodes = new_nodes
                elseif order == 3 
                    if ord_nodes[1] == nodes[1] && ord_nodes[2] == nodes[2] && ord_nodes[3] == nodes[3]
                        new_nodes = vcat(ord_nodes, nodes[4], nodes[5], nodes[8], nodes[9], nodes[6], nodes[7], nodes[10])
                    elseif ord_nodes[1] == nodes[1] && ord_nodes[2] == nodes[3] && ord_nodes[3] == nodes[2]
                        new_nodes = vcat(ord_nodes, nodes[8], nodes[9], nodes[4], nodes[5], nodes[7], nodes[6], nodes[10])
                    elseif ord_nodes[1] == nodes[2] && ord_nodes[2] == nodes[1] && ord_nodes[3] == nodes[3]
                        new_nodes = vcat(ord_nodes, nodes[5], nodes[4], nodes[6], nodes[7], nodes[8], nodes[9], nodes[10])
                    elseif ord_nodes[1] == nodes[2] && ord_nodes[2] == nodes[3] && ord_nodes[3] == nodes[1]
                        new_nodes = vcat(ord_nodes, nodes[6], nodes[7], nodes[5], nodes[4], nodes[9], nodes[8], nodes[10])
                    elseif ord_nodes[1] == nodes[3] && ord_nodes[2] == nodes[1] && ord_nodes[3] == nodes[2]
                        new_nodes = vcat(ord_nodes, nodes[9], nodes[8], nodes[7], nodes[6], nodes[4], nodes[5], nodes[10])
                    elseif ord_nodes[1] == nodes[3] && ord_nodes[2] == nodes[2] && ord_nodes[3] == nodes[1]
                        new_nodes = vcat(ord_nodes, nodes[7], nodes[6], nodes[9], nodes[8], nodes[5], nodes[4], nodes[10])
                    end
                    nodes = new_nodes
                end
                nodes = [Int64(x) for x in nodes]
                id = Int64(id)
                new(id, order, nodes, material, boundary_egdes, dof)
            end
        end

        struct My_quad <: Element
            id::Int64
            order::Int32
            nodes::Vector{Int64}
            material::Int64
            boundary_egdes::Vector{Int64}
            dof::Vector{Int64}


            # Constructor personalizado para Quad
            function My_quad(id::UInt64, order::Int32, nodes::Vector{UInt64}, material::Int64, boundary_egdes::Vector{Int64}, dof::Vector{Int64})
                if order == GMSH_QUAD_ORDER_1  
                    order = 1
                elseif order == GMSH_QUAD_ORDER_2 
                    order = 2
                else
                    throw(ArgumentError("El cuadrado debe tener orden 1 o 2"))
                end

                nodes = [Int64(x) for x in nodes]  
                id = Int64(id) 
                new(id, order, nodes, material, boundary_egdes, dof)
            end
        end

        global EPSILON_0 = 8.854e-12
        global MU_0 = 1.256e-6

        struct My_material 
            id::Int32
            descript::String
            epsilon_r::Float64
            mu_r::Float64
            epsilon_0::Float64
            mu_0::Float64

            # Constructor personalizado para My_material
            function My_material(id::Int32, descript::String, epsilon_r::Float64, mu_r::Float64)
                new(id, descript, epsilon_r, mu_r, EPSILON_0, MU_0)
            end
        end

        global DIRICHLET = 1
        global NEUMAN = 2
        global CAUCHY = 3

        struct Boundary_condition
            id::Int64
            descript::String
            type::Int64
            nodes::Vector{Int64}

            # Constructor personalizado para Boundary_condition
            function Boundary_condition(id::Int32, type::Int64, nodes::Vector{UInt64})
                if type == DIRICHLET
                    descript = "Son esenciales y deben fijarse explícitamente"
                elseif type == NEUMAN
                    descript = "Son naturales y se aplican implícitamente"
                elseif type == CAUCHY
                    descript = ""
                end
                nodes = [Int64(x) for x in nodes]
                id = convert(Int64, id)
                new(id, descript, type, nodes)
            end
        end

const to = TimerOutputs()

@timeit to "Código Completo" begin

    @timeit to "Construir y cargar la malla" begin
        @timeit to "Crear la malla desde el API de Gmsh en Julia" begin
            gmsh.initialize()

            #gmsh.model.add("rect_waveguide_triangles")

            global a = 0.02286
            global b = 0.01016
            tm = a/5

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
            gmsh.model.addPhysicalGroup(2, [1], 6, "Aire dentro de la guia")

            # Generar la malla
            gmsh.model.mesh.generate(2)

            # Generar los elementos de orden 2
            #gmsh.model.mesh.setOrder(2)

            # Obtenemos los identificadores de los nodos, las coordenadas
            nodeTags_t, coord_t = gmsh.model.mesh.getNodes()

            #= Obtenemos los tipos de elementos, sus identificadores 
            y los nodos de todos los elementos del mismo tipo concatenados,
            estos son los identificadores unicos de los nodos que forman parte de los elementos.=#
            elementTypes_t, elementTags_t, nodeTags_e_t = gmsh.model.mesh.getElements()
            element_prop_t = gmsh.model.mesh.getElementProperties(elementTypes_t[2])

            # Obtenemos las entidades
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
            groups_t = Vector{Tuple{Int32, Vector{UInt64}}}()
            for group in physical_group_t
                phy_nodes = gmsh.model.mesh.getNodesForPhysicalGroup(1, group[2])
                push!(groups_t, (group[2], phy_nodes[1]))
            end

            # Finalizar el modelo gmsh
            gmsh.finalize()
        end
        
        
        @timeit to "Construir la malla en Julia" begin
            # Vector con todas las condiciones de contorno de la malla
            mesh_cond = Vector{Boundary_condition}()

            # Vector con todos los materiales de la malla
            mesh_materials = Vector{My_material}()

            # Vector con todas las coordenadas de los nodos de la malla
            mesh_coord_2D = Vector{Float64}()

            # Vector con todos los elementos de la malla
            mesh_elements = Vector{Element}()

            @timeit to "Obtener todos los materiales" begin
                # Funcion que obtiene todos los materiales de la malla
                function get_all_materials(materials_names::Vector{Tuple{Int32, String}})
                    for materials in materials_names
                        material = My_material(materials[1], materials[2], 0.0, 0.0)
                        push!(mesh_materials, material)
                    end
                end

                get_all_materials(materials_names_t) # Clase geometria de entrada
            end
            
            @timeit to "Obtener todas las condiciones de contorno y sus nodos" begin
                # Funcion que obtiene todas las condiciones de contorno
                function get_all_boundary_cond(boundary_names::Vector{Tuple{Int32, String}}, groups::Vector{Tuple{Int32, Vector{UInt64}}})
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
                        elseif boundarys[2] == "CAUCHY"
                            type = Boundary_cond.CAUCHY
                        end
                    boundary = Boundary_cond.Boundary_condition(boundarys[1], type, nodes)
                    push!(mesh_cond, boundary)
                    end
                end

                get_all_boundary_cond(boundary_names_t, groups_t)
            end

            @timeit to "Guardas todas las coordenadas de los nodos" begin
                # Funcion que guarda las coordenadas x e y de todos los nodos de la malla en mesh_coord_2D
                function get_all_coordenates(coord::Vector{Float64})

                    # Coord contiene las coordenadas completas con {n1x, n1y, n1z, n2x, n2y, n2z, ...}

                    # Iteramos a traves de las coordenadas
                    for i in 1:3:length(coord)  # Saltamos de 3 en 3 para seleccionar x, y
                        x = coord[i]
                        y = coord[i + 1]

                        push!(mesh_coord_2D, x)
                        push!(mesh_coord_2D, y)
                    end
                end 

                get_all_coordenates(coord_t)
            end

            @timeit to "Obtener todos los elementos de la malla con sus respectivos atributos" begin
                # Funcion para determinar el numero de nodos segun el orden del elemento
                function evaluate_order(orden::Int32)
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
                function get_boundary_nodes(element::Element, boundary_conditions::Vector{Boundary_condition})
                    if isa(element, My_triangle)
                        edge_1 = [element.nodes[1], element.nodes[2]]
                        edge_2 = [element.nodes[1], element.nodes[3]]
                        edge_3 = [element.nodes[2], element.nodes[3]]
                        edges = (edge_1, edge_2, edge_3)
                    elseif isa(element, My_quad)  
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
                    for i in eachindex(elementType)
                        if elementType[i] == GMSH_TRIA_ORDER_1 || elementType[i] == GMSH_TRIA_ORDER_2 || elementType[i] == GMSH_TRIA_ORDER_3
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
                            element = My_triangle(elemtTags_e[j],elementType[i], nodesOfElement, 0, Vector{Int64}(), Vector{Int64}())

                            # Extraemos las condiciones de contorno
                            get_boundary_nodes(element, mesh_cond)

                            # Guardamos el elemento en mesh_elements
                            push!(mesh_elements, element)
                            end
                        elseif elementType[i] == GMSH_QUAD_ORDER_1 || elementType[i] == GMSH_QUAD_ORDER_2
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
                            element = My_quad(elemtTags_e[j],elementType[i], nodesOfElement, 0, Vector{Int64}(), Vector{Int64}())

                            # Extraemos las condiciones de contorno
                            get_boundary_nodes(element, mesh_cond)

                            # Guardamos el elemento en mesh_elements
                            push!(mesh_elements, element)
                            end
                        end
                    end
                end

                get_all_elements(elementTypes_t, elementTags_t, nodeTags_e_t)
            end
            
            # Funcion que guarda los dof en cada elemento
            function get_all_dofs(elements::Vector{Element}) 
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

                    if (isa(element, My_triangle) && (element.order != 1))
                        if isa(element, My_triangle)
                            edge_1 = [element.nodes[1], element.nodes[2]]
                            edge_2 = [element.nodes[1], element.nodes[3]]
                            edge_3 = [element.nodes[2], element.nodes[3]]
                            edges = (edge_1, edge_2, edge_3)

                        elseif isa(element, My_quad)
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
                    
                    if (isa(element, My_triangle) && (element.order == 3)) || (isa(element, My_quad) && (element.order == 2))
                        cont +=1
                        push!(element.dof, cont)
                    end
                end 
                return cont
            end

            total_dofs = get_all_dofs(mesh_elements)

        end
        
    end
    
    @timeit to "Calculo del método de elementos finitos con integración numérica" begin
        global FEM_TE = false

        #-------------------------------Triangulos-------------------------------#

        # Funcion que calcula las funciones de base en el triangulo de referencia en base al orden del triangulo
        # pasado por argumento
        function evalFunctionsForReference(gaussX::Float64, gaussY::Float64, element::My_triangle)

            if element.order==1
                scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY]
            elseif element.order==2
                scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY, (1-gaussX-gaussY)*gaussX, 
                (1-gaussX-gaussY)*gaussY, gaussX*gaussY]
            elseif element.order==3
                scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY,
                (1-gaussX-gaussY)*gaussX, (1-gaussX-gaussY)*gaussX*(1-2*gaussX-gaussY), 
                (1-gaussX-gaussY)*gaussY, (1-gaussX-gaussY)*gaussY*(1-gaussX-2*gaussY), 
                gaussX*gaussY , gaussX*gaussY*(gaussX-gaussY), 
                (1-gaussX-gaussY)*gaussX*gaussY]
            end

            return scalar_basis_func
        end

        # Funcion que calcula el gradiente de las funciones de base del triangulo de referencia
        function evalGradientFunction(gaussX::Float64, gaussY::Float64, element::My_triangle)

            if element.order == 1
                gradient = [-1.0 -1.0 ;
                            1.0 0.0 ;
                            0.0 1.0]
            elseif element.order == 2
                gradient = [-1.0 -1.0 ; 
                            1.0 0.0 ;
                            0.0 1.0;
                            1-2*gaussX-gaussY -gaussX ;
                            -gaussY 1-gaussX-2*gaussY ;
                            gaussY gaussX]
            elseif element.order == 3
                gradient = [-1.0 -1.0 ;
                            1.0 0.0 ;
                            0.0 1.0 ;
                            1-2*gaussX-gaussY -gaussX ;
                            6*gaussX^2-6*gaussX-2*gaussY+6*gaussX*gaussY+gaussY^2+1 3*gaussX^2+2*gaussX*gaussY-2*gaussX ;
                            -gaussY 1-gaussX-2*gaussY ;
                            3*gaussY^2+2*gaussX*gaussY-2*gaussY gaussX^2-2*gaussX-6*gaussY+6*gaussX*gaussY+6*gaussY^2+1 ;
                            gaussY gaussX ;
                            2*gaussY*gaussX-gaussY^2 gaussX^2-2*gaussX*gaussY;
                            gaussY-2*gaussY*gaussX-gaussY^2 gaussX-gaussX^2-2*gaussX*gaussY ]
            end

            return gradient
        end

        # Funcion que calcula el Jacobiano del triangulo
        function jacobian_calc(gradient::Matrix{Float64}, coord::Vector{Float64}, element::My_triangle)
            ∂x_∂ξ = 0.0
            ∂x_∂η = 0.0
            ∂y_∂ξ = 0.0
            ∂y_∂η = 0.0
            nodes = element.nodes

            rows, columns = size(gradient)

            for row in 1:rows
                for column in 1:columns
                    if (column % 2 == 0)
                        ∂x_∂η += (coord[nodes[row]*2-1] * gradient[row, column])
                        ∂y_∂η += (coord[nodes[row]*2] * gradient[row, column])
                    else
                        ∂x_∂ξ += (coord[nodes[row]*2-1] * gradient[row, column])
                        ∂y_∂ξ += (coord[nodes[row]*2] * gradient[row, column])
                    end
                end
            end

            ∂x_∂ξ_2 = -coord[nodes[1]*2-1] + coord[nodes[2]*2-1]
            ∂x_∂η_2 = -coord[nodes[1]*2-1] + coord[nodes[3]*2-1]
            ∂y_∂ξ_2 = -coord[nodes[1]*2] + coord[nodes[2]*2]
            ∂y_∂η_2 = -coord[nodes[1]*2] + coord[nodes[3]*2]

            # Construir la matriz Jacobiana
            J = [∂x_∂ξ ∂y_∂ξ; 
                ∂x_∂η  ∂y_∂η]
            J_2 = [∂x_∂ξ_2 ∂y_∂ξ_2; 
                ∂x_∂η_2  ∂y_∂η_2]
            det_J = det(J)
            det_J_2 = det(J_2)
            J_inv = inv(J)
            J_inv_2 = inv(J_2)
            
            return J, det_J_2, J_inv_2
        end

        function boundary_matrix_mass!(mass_matrix::Matrix{Float64}, element::My_triangle)

            # Comprobamos si tiene la condicion de contorno de DIRICHLET
            if (element.order == 1 || element.order == 2)
                for i in 1:length(element.boundary_egdes)
                    if element.boundary_egdes[i] == DIRICHLET
                        if element.order == 1
                            if i == 1
                                dof_index_1 = (1, 2)
                            elseif i == 2
                                dof_index_1 = (1, 3)
                            elseif i == 3
                                dof_index_1 = (2, 3)
                            end
                        elseif element.order == 2
                            if i == 1
                                dof_index_1 = (1, 2, 4)
                            elseif i == 2
                                dof_index_1 = (1, 3, 5)
                            elseif i == 3
                                dof_index_1 = (2, 3, 6)
                            end
                        end
                        
                        for dof_index in dof_index_1
                            mass_matrix[element.dof[dof_index], element.dof[dof_index]] = 1.0
                            for j in 1:total_dofs
                                if j != element.dof[dof_index]
                                    mass_matrix[element.dof[dof_index], j] = 0.0
                                    mass_matrix[j, element.dof[dof_index]] = 0.0
                                end
                            end
                        end
                    end
                end
            elseif element.order == 3
                for i in 1:length(element.boundary_egdes)
                    if element.boundary_egdes[i] == DIRICHLET
                        if i == 1
                            dof_index_1 = (1, 2, 4, 5)
                        elseif i == 2
                            dof_index_1 = (1, 3, 6, 7)
                        elseif i == 3
                            dof_index_1 = (2, 3, 7, 8)
                        end
                        
                        for dof_index in dof_index_1
                            mass_matrix[element.dof[dof_index], element.dof[dof_index]] = 1.0
                            for j in 1:total_dofs
                                if j != element.dof[dof_index]
                                    mass_matrix[element.dof[dof_index], j] = 0.0
                                    mass_matrix[j, element.dof[dof_index]] = 0.0
                                end
                            end
                        end
                    end
                end
            end
        end

        function boundary_matrix_stiffness!(stiffness_matrix::Matrix{Float64}, element::My_triangle)

            # Comprobamos si tiene la condicion de contorno de DIRICHLET
            if (element.order == 1 || element.order == 2)
                for i in 1:length(element.boundary_egdes)
                    if element.boundary_egdes[i] == DIRICHLET
                        if element.order == 1
                            if i == 1
                                dof_index_1 = (1, 2)
                            elseif i == 2
                                dof_index_1 = (1, 3)
                            elseif i == 3
                                dof_index_1 = (2, 3)
                            end
                        elseif element.order == 2
                            if i == 1
                                dof_index_1 = (1, 2, 4)
                            elseif i == 2
                                dof_index_1 = (1, 3, 5)
                            elseif i == 3
                                dof_index_1 = (2, 3, 6)
                            end
                        end
                        
                        for dof_index in dof_index_1
                            stiffness_matrix[element.dof[dof_index], element.dof[dof_index]] = 1.0e10
                            for j in 1:total_dofs
                                if j != element.dof[dof_index]
                                    stiffness_matrix[element.dof[dof_index], j] = 0.0
                                    stiffness_matrix[j, element.dof[dof_index]] = 0.0
                                end
                            end
                        end
                    end
                end
            elseif element.order == 3
                for i in 1:length(element.boundary_egdes)
                    if element.boundary_egdes[i] == DIRICHLET
                        if i == 1
                            dof_index_1 = (1, 2, 4, 5)
                        elseif i == 2
                            dof_index_1 = (1, 3, 6, 7)
                        elseif i == 3
                            dof_index_1 = (2, 3, 7, 8)
                        end
                        
                        for dof_index in dof_index_1
                            stiffness_matrix[element.dof[dof_index], element.dof[dof_index]] = 1.0e10
                            for j in 1:total_dofs
                                if j != element.dof[dof_index]
                                    stiffness_matrix[element.dof[dof_index], j] = 0.0
                                    stiffness_matrix[j, element.dof[dof_index]] = 0.0
                                end
                            end
                        end
                    end
                end
            end
        end

        #-------------------------------Matriz de masa-------------------------------#

        # Funcion que calcula la matriz de masa por elemento
        function get_local_mass_matrix(n_gauss::Matrix{Float64}, w_gauss::Vector{Float64}, element::Element)

            # Inicializar la matriz de masa
            mass_matrix = zeros(Float64, length(element.dof), length(element.dof))

            for k in 1:length(w_gauss)
                # Calculamos las funciones de base, los gradientes y el jacobiano
                basis_functions = evalFunctionsForReference(n_gauss[k,1], n_gauss[k,2], element)
                gradient = evalGradientFunction(n_gauss[k,1], n_gauss[k,2], element)
                jacobiano, det_J, inv_J = jacobian_calc(gradient, mesh_coord_2D, element)

                for i in 1:length(element.dof)
                    for j in 1:length(element.dof)
                        # Calcular la contribución al elemento de masa
                        mass_matrix[i, j] += w_gauss[k] * basis_functions[i] * basis_functions[j] * abs(det_J)        
                    end
                end
            end
            
            return mass_matrix
        end

        # Funcion que calcula la matriz de masa total
        function assemble_global_mass_matrix!(global_mass_matrix::Matrix{Float64}, local_mass_matrix::Matrix{Float64}, element::Element)
            for i in 1:length(element.dof)
                for j in 1:length(element.dof)
                    global_mass_matrix[element.dof[i], element.dof[j]] += local_mass_matrix[i, j]
                end
            end
        end

        #-------------------------------Matriz de rigidez-------------------------------#

        # Funcion que calcula la matriz de rigidez por elemento
        function get_local_stiffness_matrix(n_gauss::Matrix{Float64}, w_gauss::Vector{Float64}, element::Element)

            # Inicializar la matriz de rigidez
            stiffness_matrix = zeros(Float64, length(element.dof), length(element.dof))

            for k in 1:length(w_gauss)
                # Calculamos los gradientes y el jacobiano
                gradient = evalGradientFunction(n_gauss[k,1], n_gauss[k,2], element)
                jacobiano, det_J, inv_J = jacobian_calc(gradient, mesh_coord_2D, element)

                # Calcular la contribución al elemento de rigidez
                for i in 1:length(element.dof)
                    ∇Φ_i = [gradient[i, 1] gradient[i, 2]]  # Obtiene el gradiente de la i-ésima función de forma
                    for j in 1:length(element.dof)
                        ∇Φ_j = [gradient[j, 1] gradient[j, 2]]  # Obtiene el gradiente de la j-ésima función de forma
                        stiffness_matrix[i,j] += w_gauss[k] * (inv_J * ∇Φ_i') ⋅ (inv_J * ∇Φ_j') * abs(det_J)
                    end
                end        
            end

            return stiffness_matrix
        end

        # Funcion que calcula la matriz de rigidez total
        function assemble_global_stiffness_matrix!(global_stiffness_matrix::Matrix{Float64}, local_stiffness_matrix::Matrix{Float64}, element::Element)

            # Ensamblar la matriz de rigidez local en la matriz de rigidez global
            for i in 1:length(element.dof)
                for j in 1:length(element.dof)
                    global_stiffness_matrix[element.dof[i], element.dof[j]] += local_stiffness_matrix[i, j]
                end
            end

            return global_stiffness_matrix
        end

        #-------------------------------Calculos-------------------------------#

        # Definimos las matrices globales y las calculamos
        global_mass_matrix = zeros(Float64, total_dofs, total_dofs)
        global_stiffness_matrix = zeros(Float64, total_dofs, total_dofs)

        # Vertices del triangulo de referencia
        vertex = [
            0.0 0.0
            1.0 0.0
            0.0 1.0
        ]

        # Creamos las matrices de rigidez y masa globales
        for element in mesh_elements

            # Obtener los puntos y los pesos de la regla de cuadratura Gauss-Legendre
            n_gauss, w_gauss = simplexquad(6, vertex)

            @timeit "Calculo de la matriz de masa global" begin
                local_mass_matrix = get_local_mass_matrix(n_gauss, w_gauss, element)
                assemble_global_mass_matrix!(global_mass_matrix, local_mass_matrix, element)
            end
            
            
            @timeit "Calculo de la matriz de rigidez global" begin
                local_stiffness_matrix = get_local_stiffness_matrix(n_gauss, w_gauss, element)
                assemble_global_stiffness_matrix!(global_stiffness_matrix, local_stiffness_matrix, element)
            end
            
            
            
        end

        # Aplicamos las condiciones de contorno de DIRICHLET
        @timeit "Aplicación de las condiciones de contorno de DIRICHLET" begin
            if FEM_TE == false
                for element in mesh_elements
                    boundary_matrix_mass!(global_mass_matrix, element)
                
                    boundary_matrix_stiffness!(global_stiffness_matrix, element) 
                end
            end
        end
        

        @timeit "Calculo del error relativo" begin
            eigenv = eigen( global_stiffness_matrix, global_mass_matrix)

            k_c_fem = real(eigenv.values)

            # Definimos el k_c^2 analitico
            k_c_analytic = Vector{Float64}()

            # Calculamos los k_C^2 analiticos
            k_c_analytic_01 = (((0*pi)/a)^2+((1*pi)/b)^2)
            k_c_analytic_10 = (((1*pi)/a)^2+((0*pi)/b)^2)
            k_c_analytic_11 = (((1*pi)/a)^2+((1*pi)/b)^2)
            push!(k_c_analytic, k_c_analytic_01)
            push!(k_c_analytic, k_c_analytic_10)
            push!(k_c_analytic, k_c_analytic_11)
            sort!(k_c_analytic)

            # Funcion que calcula el error relativo
            rel_error = Vector{Float64}()
            function rel_err!(k_c_fem::Vector{Float64}, k_c_analytic::Vector{Float64}, rel_error::Vector{Float64})

                if FEM_TE == true
                    error = abs((k_c_fem[2]-k_c_analytic[1])/k_c_analytic[1])
                else
                    error = abs((k_c_fem[1]-k_c_analytic[3])/k_c_analytic[3])
                end
                
                push!(rel_error, error)
                
            end

            rel_err!(k_c_fem, k_c_analytic, rel_error)

        end
    end
end





