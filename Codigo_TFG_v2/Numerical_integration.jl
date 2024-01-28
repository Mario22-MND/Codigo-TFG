module Numerical_integration

using LinearAlgebra

include("Mesh_object.jl")

#-------------------------------Triangulos-------------------------------#

# Funcion que calcula las funciones de base en el triangulo de referencia en base al orden del triangulo
# pasado por argumento
function evalFunctionsForReferenceTriangle(gaussX::Float64, gaussY::Float64, order::Int32)
    scalar_basis_func = Vector{Float64}()

    if element.order==1
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY]
    elseif element.order==2
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY, (1-gaussX-gaussY)*gaussX, 
        (1-gaussX-gaussY)*gaussY, gaussX*gaussY]
    elseif element.order==3
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY, (1-gaussX-gaussY)*gaussX*(1-2*gaussX-gaussY), 
        (1-gaussX-gaussY)*gaussY*(1-gaussX-2*gaussY), gaussX*gaussY*(gaussX-gaussY), (1-gaussX-gaussY)*gaussX*(1-2*gaussX-gaussY), 
        (1-gaussX-gaussY)*gaussY*(1-gaussX-2*gaussY), gaussX*gaussY*(gaussX-gaussY),
        (1-gaussX-gaussY)*gaussX*gaussY]
    end
    
    return scalar_basis_func
end

# Funcion que calcula el gradiente de las funciones de base del triangulo de referencia
function evalGradientTriangleFunction(gaussX::Float64, gaussY::Float64, order::Int32, scalar_basis_func::Vector{Float64})
    gradient = Vector{Float64}(undef, length(scalar_basis_func)*2)

    if element.order==1
        gradient = [-1, -1, 1, 0, 0, 1]
    elseif element.order==2
        gradient = [-1, -1, 1, 0, 0, 1, 1-2*gaussX-gaussY, -gaussX, -gaussY, 1-gaussX-2*gaussY, 
        gaussY, gaussX]
    elseif element.order==3
        gradient = [-1.0, -1.0, 1.0, 0.0, 0.0, 1.0,
        6*gaussX^2-6*gaussX-2*gaussY+6*gaussX*gaussY+gaussY^2+1, 3*gaussX^2-2*gaussX*gaussY-2*gaussX, 
        6*gaussX^2-6*gaussX-2*gaussY+6*gaussX*gaussY+gaussY^2+1, 3*gaussX^2-2*gaussX*gaussY-2*gaussX,
        3*gaussY^2-2*gaussX*gaussY-2*gaussY, gaussX^2-2*gaussX-6*gaussY+6*gaussX*gaussY+6*gaussY^2+1, 
        3*gaussY^2-2*gaussX*gaussY-2*gaussY, gaussX^2-2*gaussX-6*gaussY+6*gaussX*gaussY+6*gaussY^2+1,
        2*gaussY*gaussX-gaussY^2, gaussX^2-2*gaussX*gaussY, 
        2*gaussY*gaussX-gaussY^2, gaussX^2-2*gaussX*gaussY, 
        gaussY-2*gaussY*gaussX-gaussY^2, gaussX-gaussX^2-2*gaussX*gaussY]
    end

    return gradient
end

# Funcion que calcula el Jacobiano del triangulo
function jacobian_calc_t(gradient::Vector{Float64}, coord::Vector{Float64}, nodes::Vector{Int64})
    ∂x_∂ξ = 0.0
    ∂x_∂η = 0.0
    ∂y_∂ξ = 0.0
    ∂y_∂η = 0.0
    if size(gradient) == 14

        vertex = nodes[1:3]
        edges_1 = nodes[4:2:9]
        edges_2 = nodes[5:2:10]
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

#-------------------------------Cuadrilateros-------------------------------#

# Funcion que calcula las funciones de base en el cuadrilatero de referencia en base al orden del 
# cuadrilatero pasado por argumento
function evalFunctionsForReferenceQuad(gaussX::Float64, gaussY::Float64, order::Int32)
    scalar_basis_func = Vector{Float64}()
    if order==1
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY]
    elseif order==2
        scalar_basis_func = [1-gaussX-gaussY, gaussX, gaussY, (1-gaussX-gaussY)*gaussX, 
        (1-gaussX-gaussY)*gaussY, gaussX*gaussY]
    end
    
    return scalar_basis_func
end

# Funcion que calcula el gradiente de las funciones de base del cuadrilatero de referencia
function evalGradientQuadFunction(gaussX::Float64, gaussY::Float64, order::Int32, scalar_basis_func::Vector{Float64})
    gradient = Vector{Float64}(undef, length(scalar_basis_func)*2)

    if order==1
        gradient = [-1, -1, 1, 0, 0, 1]
    elseif order==2
        gradient = [-1, -1, 1, 0, 0, 1, 1-2*gaussX-gaussY, -gaussX, -gaussY, 1-gaussX-2*gaussY, 
        gaussY, gaussX]
    end

    return gradient
end

# Funcion que calcula el Jacobiano del cuadrilatero
function jacobian_calc_c(gradient::Vector{Float64}, coord::Vector{Float64}, nodes::Vector{Int64})
    ∂x_∂ξ = 0.0
    ∂x_∂η = 0.0
    ∂y_∂ξ = 0.0
    ∂y_∂η = 0.0
    if size(gradient) == 14

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

#-------------------------------Calculos-------------------------------#

gaussX = 1.0
gaussY = 1.0

for element in Mesh_object.mesh_elements
    if isa(element, Mesh_object.Elements.My_triangle)
        basis_function = evalFunctionsForReferenceTriangle( gaussX, gaussY, element.order)
        gradient = evalGradientTriangleFunction(gaussX, gaussY, element.order, basis_function)
        jacobiano, det_J, inv_J = jacobian_calc_t(gradient, Mesh_object.mesh_coord_2D, element.nodes)
    elseif isa(element, Mesh_object.Elements.My_quad)
        basis_function = evalFunctionsForReferenceQuad( gaussX, gaussY, element.order)
        gradient = evalGradientQuadFunction(gaussX, gaussY, element.order, basis_function)
        jacobiano, det_J, inv_J = jacobian_calc_c(gradient, Mesh_object.mesh_coord_2D, element.nodes)
    end

end

end