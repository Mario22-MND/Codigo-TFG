using SimplexQuad

# Definir la función a integrar
#f(x, y) = 1-x-y*x
f(x, y) = x

vertices = [
    0.0 0.0
    1.0 0.0
    0.0 1.0
]

npoints = 2

X, W = simplexquad(npoints, vertices)

result = sum(W[i] * f(X[i,1],X[i,2]) for i in 1:length(W))

# Función para crear el modo con el formato (a, b)
function crear_modo(a, b)
    return "Modo ($a,$b)"
end

# Vector para almacenar los modos
modos = []

# Bucle para generar y almacenar los modos
for i in 0:2
    for j in 0:2
        push!(modos, crear_modo(i, j))
    end
end

# Mostrar los modos generados
for modo in modos
    println(modo)
end

matriz = [1 2 ; 3 4]

abs(-1)