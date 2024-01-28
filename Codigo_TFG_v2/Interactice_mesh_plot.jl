include("Mesh_object.jl")

using Plots
using WebIO

coord = Mesh_object.mesh_coord_2D
elements = Mesh_object.mesh_elements
b_conditions = Mesh_object.mesh_cond
order = 0

# Variable para cambiar de mallas triangulares a mallas cuadradas
quad_meshes = Mesh_object.quadrangular_meshes;

# Extraemos todas las coordenadas x e y separadas
coord_x = [coord[i] for i in 1:2:length(coord)]
coord_y = [coord[i] for i in 2:2:length(coord)]

# Encontramos los valores minimos y maximos de x e y
min_x = minimum(coord_x)
max_x = maximum(coord_x)
min_y = minimum(coord_y)
max_y = maximum(coord_y)

# Creamos un rango de valores para los ejes x e y con una separacion de 1 punto
x_range = min_x:1:max_x
y_range = min_y:1:max_y


# Creamos un grafico de dispersion
# Crear un gráfico interactivo
gr()
scatter(coord_x, coord_y, legend=false, xlabel="X", ylabel="Y", xticks=x_range, yticks=y_range)

# Iteramos a traves de los elementos para pintarlos
for element in elements
    element_nodes = element.nodes
    order = element.order

    # Filtrar solo los nodos que son los vertices del triangulo o cuadrilatero
    if quad_meshes == false
        # Dibujamos los triangulos
        vertex_nodes = element_nodes[1:3]
    else
        # Dibujamos los cuadrilateros
        vertex_nodes = element_nodes[1:4]
    end
    

    # Dibuja un triangulo utilizando solo los nodos de los vertices
    x = [coord[node * 2 - 1] for node in vertex_nodes]
    y = [coord[node * 2] for node in vertex_nodes]

    # Dibuja los elementos
    plot!(x, y, seriestype = :shape, linecolor = :blue, fillalpha = 0)
end

# Iteramos a traves de las aristas con condicion de contorno para pintarlas
for element in elements
    nodes = element.nodes
    b_edges = element.boundary_egdes
    if quad_meshes == false
        if b_edges[1] != 0
            x = [coord[nodes[1] * 2 - 1], coord[nodes[2] * 2 - 1]]
            y = [coord[nodes[1] * 2], coord[nodes[2] * 2]]
            for b_cond in b_conditions
                if b_edges[1] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
        if b_edges[2] != 0
            x = [coord[nodes[1] * 2 - 1], coord[nodes[3] * 2 - 1]]
            y = [coord[nodes[1] * 2], coord[nodes[3] * 2]]
            for b_cond in b_conditions
                if b_edges[2] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
        if b_edges[3] != 0
            x = [coord[nodes[2] * 2 - 1], coord[nodes[3] * 2 - 1]]
            y = [coord[nodes[2] * 2], coord[nodes[3] * 2]]
            for b_cond in b_conditions
                if b_edges[3] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
    else
        if b_edges[1] != 0
            x = [coord[nodes[1] * 2 - 1], coord[nodes[2] * 2 - 1]]
            y = [coord[nodes[1] * 2], coord[nodes[2] * 2]]
            for b_cond in b_conditions
                if b_edges[1] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
        if b_edges[2] != 0
            x = [coord[nodes[2] * 2 - 1], coord[nodes[3] * 2 - 1]]
            y = [coord[nodes[2] * 2], coord[nodes[3] * 2]]
            for b_cond in b_conditions
                if b_edges[2] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
        if b_edges[3] != 0
            x = [coord[nodes[3] * 2 - 1], coord[nodes[4] * 2 - 1]]
            y = [coord[nodes[3] * 2], coord[nodes[4] * 2]]
            for b_cond in b_conditions
                if b_edges[3] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
        if b_edges[4] != 0
            x = [coord[nodes[1] * 2 - 1], coord[nodes[4] * 2 - 1]]
            y = [coord[nodes[1] * 2], coord[nodes[4] * 2]]
            for b_cond in b_conditions
                if b_edges[4] == b_cond.id
                    if b_cond.type == Mesh_object.Boundary_cond.DIRICHLET
                        plot!(x, y, line=:solid, linecolor = :red, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.NEUMAN
                        plot!(x, y, line=:solid, linecolor = :black, fillalpha = 0, linewidth=3)
                    elseif b_cond.type == Mesh_object.Boundary_cond.CAUCHY
                        plot!(x, y, line=:solid, linecolor = :green, fillalpha = 0, linewidth=3)
                    end
                end
            end
        end
    end
end



if quad_meshes == false 
    if order == 1 
        if length(b_conditions) == 1
            savefig("Mesh_interactive\\triangular_mesh_order_1_boundary_1.html")
            display(PlotlyJS.open("output_plot.html"))
        elseif length(b_conditions) == 2
            savefig("Mesh_interactive\\triangular_mesh_order_1_boundary_2.html")
            display(PlotlyJS.open("output_plot.html"))
        elseif length(b_conditions) == 3
            savefig("Mesh_interactive\\triangular_mesh_order_1_boundary_3.html")
            WebIO.open("Mesh_interactive\\triangular_mesh_order_1_boundary_3.html")
        end
    elseif order == 2
        if length(b_conditions) == 1
            savefig("Mesh_interactive\\triangular_mesh_order_2_boundary_1.html")
            display(PlotlyJS.open("output_plot.html"))
        elseif length(b_conditions) == 2
            savefig("Mesh_interactive\\triangular_mesh_order_2_boundary_2.html")
            display(PlotlyJS.open("output_plot.html"))
        elseif length(b_conditions) == 3
            savefig("Mesh_interactive\\triangular_mesh_order_2_boundary_3.html")
            display(PlotlyJS.open("output_plot.html"))
        end
    elseif order == 3
        if length(b_conditions) == 1
            savefig("Mesh_interactive\\triangular_mesh_order_3_boundary_1.html")
            display(PlotlyJS.open("output_plot.html"))
        elseif length(b_conditions) == 2
            savefig("Mesh_interactive\\triangular_mesh_order_3_boundary_2.html")
            display(PlotlyJS.open("output_plot.html"))
        elseif length(b_conditions) == 3
            savefig("Mesh_interactive\\triangular_mesh_order_3_boundary_3.html")
            display(PlotlyJS.open("output_plot.html"))
        end
    end
else
    if order == 1 
        if length(b_conditions) == 1
            savefig("Mesh_interactive\\quad_mesh_order_1_boundary_1.html")
            display(PlotlyJS.open("Mesh_interactive\\quad_mesh_order_1_boundary_1.html"))
        elseif length(b_conditions) == 2
            savefig("Mesh_interactive\\quad_mesh_order_1_boundary_2.html")
            display(PlotlyJS.open("Mesh_interactive\\quad_mesh_order_1_boundary_2.html"))
        elseif length(b_conditions) == 3
            savefig("Mesh_interactive\\quad_mesh_order_1_boundary_3.html")
            display(PlotlyJS.open("Mesh_interactive\\quad_mesh_order_1_boundary_3.html"))
        end
    elseif order == 2
        if length(b_conditions) == 1
            savefig("Mesh_interactive\\quad_mesh_order_2_boundary_1.html")
            display(PlotlyJS.open("Mesh_interactive\\quad_mesh_order_2_boundary_1.html"))
        elseif length(b_conditions) == 2
            savefig("Mesh_interactive\\quad_mesh_order_2_boundary_2.html")
            display(PlotlyJS.open("Mesh_interactive\\quad_mesh_order_2_boundary_2.html"))
        elseif length(b_conditions) == 3
            savefig("Mesh_interactive\\quad_mesh_order_2_boundary_3.html")
            display(PlotlyJS.open("Mesh_interactive\\quad_mesh_order_2_boundary_3.html"))
        end
    end
end

