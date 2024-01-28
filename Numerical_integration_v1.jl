module Numerical_integration

using LinearAlgebra
using SimplexQuad

include("Mesh_object.jl")

global FEM_TE = false

#-------------------------------Triangulos-------------------------------#

# Funcion que calcula las funciones de base en el triangulo de referencia en base al orden del triangulo
# pasado por argumento
function evalFunctionsForReference(gaussX::Float64, gaussY::Float64, element::Mesh_object.Elements.My_triangle)

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
function evalGradientFunction(gaussX::Float64, gaussY::Float64, element::Mesh_object.Elements.My_triangle)

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
function jacobian_calc(gradient::Matrix{Float64}, coord::Vector{Float64}, element::Mesh_object.Elements.My_triangle)
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

    #=∂x_∂ξ_2 = -coord[nodes[1]*2-1] + coord[nodes[2]*2-1] + coord[nodes[4]*2-1] * gradient[4,1] + coord[nodes[5]*2-1] * gradient[5,1] + coord[nodes[6]*2-1] * gradient[6,1]
    ∂x_∂η_2 = -coord[nodes[1]*2-1] + coord[nodes[3]*2-1] + coord[nodes[4]*2-1] * gradient[4,2] + coord[nodes[5]*2-1] * gradient[5,2] +  coord[nodes[6]*2-1] * gradient[6,2]
    ∂y_∂ξ_2 = -coord[nodes[1]*2] + coord[nodes[2]*2] + coord[nodes[4]*2] * gradient[4,1] + coord[nodes[5]*2] * gradient[5,1] + coord[nodes[6]*2] * gradient[6,1]
    ∂y_∂η_2 = -coord[nodes[1]*2] + coord[nodes[3]*2] + coord[nodes[4]*2] * gradient[4,2] + coord[nodes[5]*2] * gradient[5,2] + coord[nodes[6]*2] * gradient[6,2]=#

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

function boundary_matrix_mass!(mass_matrix::Matrix{Float64}, element::Mesh_object.Elements.My_triangle)

   # Comprobamos si tiene la condicion de contorno de DIRICHLET
    if (element.order == 1 || element.order == 2)
        for i in 1:length(element.boundary_egdes)
            if element.boundary_egdes[i] == Mesh_object.Boundary_cond.DIRICHLET
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
                    for j in 1:Mesh_object.total_dofs
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
            if element.boundary_egdes[i] == Mesh_object.Boundary_cond.DIRICHLET
                if i == 1
                    dof_index_1 = (1, 2, 4, 5)
                elseif i == 2
                    dof_index_1 = (1, 3, 6, 7)
                elseif i == 3
                    dof_index_1 = (2, 3, 7, 8)
                end
                
                for dof_index in dof_index_1
                    mass_matrix[element.dof[dof_index], element.dof[dof_index]] = 1.0
                    for j in 1:Mesh_object.total_dofs
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

function boundary_matrix_stiffness!(stiffness_matrix::Matrix{Float64}, element::Mesh_object.Elements.My_triangle)

    # Comprobamos si tiene la condicion de contorno de DIRICHLET
    if (element.order == 1 || element.order == 2)
        for i in 1:length(element.boundary_egdes)
            if element.boundary_egdes[i] == Mesh_object.Boundary_cond.DIRICHLET
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
                    for j in 1:Mesh_object.total_dofs
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
            if element.boundary_egdes[i] == Mesh_object.Boundary_cond.DIRICHLET
                if i == 1
                    dof_index_1 = (1, 2, 4, 5)
                elseif i == 2
                    dof_index_1 = (1, 3, 6, 7)
                elseif i == 3
                    dof_index_1 = (2, 3, 7, 8)
                end
                
                for dof_index in dof_index_1
                    stiffness_matrix[element.dof[dof_index], element.dof[dof_index]] = 1.0e10
                    for j in 1:Mesh_object.total_dofs
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

#-------------------------------Cuadrilateros-------------------------------#

# Funcion que calcula las funciones de base en el cuadrilatero de referencia en base al orden del 
# cuadrilatero pasado por argumento
function evalFunctionsForReference(gaussX::Float64, gaussY::Float64, element::Mesh_object.Elements.My_quad)
    scalar_basis_func = Vector{Float64}()
    if element.order==1
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY]
    elseif element.order==2
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY, (1-gaussX-gaussY)*gaussX, 
        (1-gaussX-gaussY)*gaussY, gaussX*gaussY]
    end
    return scalar_basis_func
end

# Funcion que calcula el gradiente de las funciones de base del cuadrilatero de referencia
function evalGradientFunction(gaussX::Float64, gaussY::Float64, scalar_basis_func::Vector{Float64}, element::Mesh_object.Elements.My_quad)
    gradient = Vector{Float64}(undef, length(scalar_basis_func)*2)

    if element.order==1
        gradient = [-1, -1, 1, 0, 0, 1]
    elseif element.order==2
        gradient = [-1, -1, 1, 0, 0, 1, 1-2*gaussX-gaussY, -gaussX, -gaussY, 1-gaussX-2*gaussY, 
        gaussY, gaussX]
    end
    return gradient
end

# Funcion que calcula el Jacobiano del cuadrilatero
function jacobian_calc(gradient::Vector{Float64}, coord::Vector{Float64}, element::Mesh_object.Elements.My_quad)
    ∂x_∂ξ = 0.0
    ∂x_∂η = 0.0
    ∂y_∂ξ = 0.0
    ∂y_∂η = 0.0
    nodes = element.nodes
    if element.order == 2

        vertex = nodes[1:3]
        edges_1 = nodes[4:2:9]
        edges_2 = nodes[4:1:9]
        face = nodes[10]

        # Calcular las derivadas parciales de x e y con respecto a xi y eta en los vertices
        for i in 1:2:length(gradient[1:6])
            ∂x_∂ξ += (coord[vertex[Int((i+1)/2)]*2-1] * gradient[i])
            ∂y_∂ξ += (coord[vertex[Int((i+1)/2)]*2] * gradient[i])
        end

        for i in 2:2:length(gradient[1:6])
            ∂x_∂η += (coord[vertex[Int(i/2)]*2-1] * gradient[i])
            ∂y_∂η += (coord[vertex[Int(i/2)]*2] * gradient[i])
        end

        # Calcular las derivadas parciales de x e y con respecto a xi y eta en las aristas
        for i in 1:2:length(gradient[7:12])
            ∂x_∂ξ += (coord[edges_1[Int((i+1)/2)]*2-1] * gradient[i]) + (coord[edges_2[Int((i+1)/2)]*2-1] * gradient[i])
            ∂y_∂ξ += (coord[edges_1[Int((i+1)/2)]*2] * gradient[i]) + (coord[edges_2[Int((i+1)/2)]*2] * gradient[i])
        end

        for i in 2:2:length(gradient[7:12])
            ∂x_∂η += (coord[edges_1[Int(i/2)]*2-1] * gradient[i]) + (coord[edges_2[Int(i/2)]*2-1] * gradient[i])
            ∂y_∂η += (coord[edges_1[Int(i/2)]*2] * gradient[i]) + (coord[edges_2[Int(i/2)]*2] * gradient[i])
        end

        # Calcular las derivadas parciales de x e y con respecto a xi y eta en la cara
        ∂x_∂ξ += (coord[face*2-1] * gradient[13])
        ∂y_∂ξ += (coord[face*2] * gradient[13])
        
        ∂x_∂η += (coord[face*2-1] * gradient[14])
        ∂y_∂η += (coord[face*2] * gradient[14])
    else

        # Calcular las derivadas parciales de x e y con respecto a xi y eta
        for i in 1:2:length(gradient)
            ∂x_∂ξ += (coord[nodes[Int((i+1)/2)]*2-1] * gradient[i])
            ∂y_∂ξ += (coord[nodes[Int((i+1)/2)]*2] * gradient[i])
        end

        for i in 2:2:length(gradient)
            ∂x_∂η += (coord[nodes[Int(i/2)]*2-1] * gradient[i])
            ∂y_∂η += (coord[nodes[Int(i/2)]*2] * gradient[i])
        end
    end

    # Construir la matriz Jacobiana
    J = [∂x_∂ξ ∂x_∂η; ∂y_∂ξ ∂y_∂η]
    det_J = det(J)
    J_inv = inv(J)
    
    return J, det_J, J_inv
end

#-------------------------------Matriz de masa-------------------------------#

# Funcion que calcula la matriz de masa por elemento
function get_local_mass_matrix(n_gauss::Matrix{Float64}, w_gauss::Vector{Float64}, element::Mesh_object.Elements.Element)

    # Inicializar la matriz de masa
    mass_matrix = zeros(Float64, length(element.dof), length(element.dof))

    for k in 1:length(w_gauss)
        # Calculamos las funciones de base, los gradientes y el jacobiano
        basis_functions = evalFunctionsForReference(n_gauss[k,1], n_gauss[k,2], element)
        gradient = evalGradientFunction(n_gauss[k,1], n_gauss[k,2], element)
        jacobiano, det_J, inv_J = jacobian_calc(gradient, Mesh_object.mesh_coord_2D, element)

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
function assemble_global_mass_matrix!(global_mass_matrix::Matrix{Float64}, local_mass_matrix::Matrix{Float64}, element::Mesh_object.Elements.Element)
    for i in 1:length(element.dof)
        for j in 1:length(element.dof)
            global_mass_matrix[element.dof[i], element.dof[j]] += local_mass_matrix[i, j]
        end
    end
end

#-------------------------------Matriz de rigidez-------------------------------#

# Funcion que calcula la matriz de rigidez por elemento
function get_local_stiffness_matrix(n_gauss::Matrix{Float64}, w_gauss::Vector{Float64}, element::Mesh_object.Elements.Element)

    # Inicializar la matriz de rigidez
    stiffness_matrix = zeros(Float64, length(element.dof), length(element.dof))

    for k in 1:length(w_gauss)
        # Calculamos los gradientes y el jacobiano
        gradient = evalGradientFunction(n_gauss[k,1], n_gauss[k,2], element)
        jacobiano, det_J, inv_J = jacobian_calc(gradient, Mesh_object.mesh_coord_2D, element)

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
function assemble_global_stiffness_matrix!(global_stiffness_matrix::Matrix{Float64}, local_stiffness_matrix::Matrix{Float64}, element::Mesh_object.Elements.Element)

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
global_mass_matrix = zeros(Float64, Mesh_object.total_dofs, Mesh_object.total_dofs)
global_stiffness_matrix = zeros(Float64, Mesh_object.total_dofs, Mesh_object.total_dofs)

# Vertices del triangulo de referencia
vertex = [
    0.0 0.0
    1.0 0.0
    0.0 1.0
]

# Creamos las matrices de rigidez y masa globales
for element in Mesh_object.mesh_elements

    # Obtener los puntos y los pesos de la regla de cuadratura Gauss-Legendre
    n_gauss, w_gauss = simplexquad(6, vertex)

    local_mass_matrix = get_local_mass_matrix(n_gauss, w_gauss, element)
    assemble_global_mass_matrix!(global_mass_matrix, local_mass_matrix, element)
    

    local_stiffness_matrix = get_local_stiffness_matrix(n_gauss, w_gauss, element)
    assemble_global_stiffness_matrix!(global_stiffness_matrix, local_stiffness_matrix, element)
    
    
end

# Aplicamos las condiciones de contorno de DIRICHLET
if FEM_TE == false
    for element in Mesh_object.mesh_elements
        boundary_matrix_mass!(global_mass_matrix, element)
    
        boundary_matrix_stiffness!(global_stiffness_matrix, element) 
    end
end


eigenv = eigen( global_stiffness_matrix, global_mass_matrix)

# Para comprobar
val_mass = real(eigvals(global_mass_matrix))
val_stiff = real(eigvals(global_stiffness_matrix))

k_c_fem = real(eigenv.values)

# Definimos el k_c^2 analitico
k_c_analytic = Vector{Float64}()

# Funcion que calcula los k_c^2 analiticos hasta los modos 4,4
function k_analytic(k_c_analytic::Vector{Float64})
    for i in 0:3
        for j in 0:3
            k_c_analytic_i = (((i*pi)/Mesh_object.Rect_waveguide_triangles.a)^2+((j*pi)/Mesh_object.Rect_waveguide_triangles.b)^2)
            push!(k_c_analytic, k_c_analytic_i)
        end
    end

    length_kc_analytic = length(k_c_analytic)
    k_c_analytic = k_c_analytic[2:length_kc_analytic]
    #= k_c_analytic = [k_c_analytic[1], k_c_analytic[5], k_c_analytic[6], k_c_analytic[6],
                    k_c_analytic[2], k_c_analytic[10], k_c_analytic[7], k_c_analytic[11], k_c_analytic[12], k_c_analytic[12],
                    k_c_analytic[3], k_c_analytic[15], k_c_analytic[8], k_c_analytic[16], k_c_analytic[13], k_c_analytic[17], k_c_analytic[18], k_c_analytic[18],
                    k_c_analytic[4], k_c_analytic[20], k_c_analytic[9], k_c_analytic[21], k_c_analytic[14], k_c_analytic[22], k_c_analytic[19], k_c_analytic[23], k_c_analytic[24], k_c_analytic[24]]
    =#
    k_c_analytic = sort(k_c_analytic)
    return k_c_analytic
end



# Calculamos los k_C^2 analiticos
#k_c_analytic = k_analytic(k_c_analytic)
k_c_analytic_01 = (((0*pi)/Mesh_object.Rect_waveguide_triangles.a)^2+((1*pi)/Mesh_object.Rect_waveguide_triangles.b)^2)
k_c_analytic_10 = (((1*pi)/Mesh_object.Rect_waveguide_triangles.a)^2+((0*pi)/Mesh_object.Rect_waveguide_triangles.b)^2)
k_c_analytic_11 = (((1*pi)/Mesh_object.Rect_waveguide_triangles.a)^2+((1*pi)/Mesh_object.Rect_waveguide_triangles.b)^2)
push!(k_c_analytic, k_c_analytic_01)
push!(k_c_analytic, k_c_analytic_10)
push!(k_c_analytic, k_c_analytic_11)
sort!(k_c_analytic)

# Funcion que calcula el error relativo
rel_error = Vector{Float64}()
function rel_err!(k_c_fem::Vector{Float64}, k_c_analytic::Vector{Float64}, rel_error::Vector{Float64})

    #=l_analytic = length(k_c_analytic)
    l_fem = length(k_c_fem)

    l_total = 0
    if l_fem >= l_analytic
        l_total = l_analytic
    else
        l_total = l_fem
    end

    for i in 1:l_total
        error = (k_c_fem[i]-k_c_analytic[i])/k_c_analytic[i]
        push!(rel_error, error)
    end=#
    if FEM_TE == true
        error = abs((k_c_fem[2]-k_c_analytic[1])/k_c_analytic[1])
    else
        error = abs((k_c_fem[1]-k_c_analytic[3])/k_c_analytic[3])
    end
    
    push!(rel_error, error)
    
end

rel_err!(k_c_fem, k_c_analytic, rel_error)

#= Comprobaciones

b = ones(Float64, Mesh_object.total_dofs)

x = inv(global_mass_matrix)*b
y = inv(global_stiffness_matrix)*b =#

end