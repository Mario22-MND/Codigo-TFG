module Elements

include("Elements_order.jl")

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
        if order == Elements_order.GMSH_TRIA_ORDER_1
            order = 1
        elseif order == Elements_order.GMSH_TRIA_ORDER_2 
            order = 2
        elseif order == Elements_order.GMSH_TRIA_ORDER_3
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
        if order == Elements_order.GMSH_QUAD_ORDER_1  
            order = 1
        elseif order == Elements_order.GMSH_QUAD_ORDER_2 
            order = 2
        else
            throw(ArgumentError("El cuadrado debe tener orden 1 o 2"))
        end

        #= Determinar las aristas de los elementos
        ord_nodes = sort(nodes[1:4])
        if order == 1
            nodes = ord_nodes
        elseif order == 2
            if ord_nodes[1] == nodes[1] 
                if ord_nodes[2] == nodes[4]
                    new_nodes = vcat(ord_nodes, nodes[6], nodes[7], nodes[8], nodes[5], nodes[9])
                elseif ord_nodes[2] == nodes[2]
                    new_nodes = vcat(ord_nodes, nodes[5], nodes[8], nodes[7], nodes[6], nodes[9])
                end
            elseif ord_nodes[1] == nodes[4]  
                if ord_nodes[2] == nodes[3]
                    new_nodes = vcat(ord_nodes, nodes[7], nodes[8], nodes[5], nodes[6], nodes[9])
                elseif ord_nodes[2] == nodes[1]
                    new_nodes = vcat(ord_nodes, nodes[6], nodes[5], nodes[8], nodes[7], nodes[9])
                end
            elseif ord_nodes[1] == nodes[3]  
                if ord_nodes[2] == nodes[2]
                    new_nodes = vcat(ord_nodes, nodes[8], nodes[5], nodes[6], nodes[7], nodes[9])
                elseif ord_nodes[2] == nodes[4]
                    new_nodes = vcat(ord_nodes, nodes[7], nodes[6], nodes[5], nodes[8], nodes[9])
                end
            elseif ord_nodes[1] == nodes[2] 
                if ord_nodes[2] == nodes[1]
                    new_nodes = vcat(ord_nodes, nodes[5], nodes[6], nodes[7], nodes[8], nodes[9])
                elseif ord_nodes[2] == nodes[3]
                    new_nodes = vcat(ord_nodes, nodes[8], nodes[7], nodes[6], nodes[5], nodes[9])
                end
            end
            nodes = new_nodes
        end=#
        nodes = [Int64(x) for x in nodes]  
        id = Int64(id) 
        new(id, order, nodes, material, boundary_egdes, dof)
    end
end

end