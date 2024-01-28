using Plots

include("Numerical_integration_v1.jl")

# Todo con orden de integracon 6

# Cargamos los resultados de orden 1
k_c_fem_TE_1 = Vector{Float64}() 
push!(k_c_fem_TE_1, 19307.37461259824)
push!(k_c_fem_TE_1, 19033.165264813615)
push!(k_c_fem_TE_1, 18996.141915532568)
push!(k_c_fem_TE_1, 18915.09942826731)

k_c_fem_TM_1 = Vector{Float64}()
push!(k_c_fem_TM_1, 126193.35118204288)
push!(k_c_fem_TM_1, 121097.32528778014)
push!(k_c_fem_TM_1, 118560.57401613961)
push!(k_c_fem_TM_1, 115615.29430584713)


# Cargamos los resultados de orden 2
k_c_fem_TE_2 = Vector{Float64}()
push!(k_c_fem_TE_2, 18889.332848919395)
push!(k_c_fem_TE_2, 18886.606021718842)
push!(k_c_fem_TE_2, 18886.445895160115)
push!(k_c_fem_TE_2, 18886.326007928707)

k_c_fem_TM_2 = Vector{Float64}() 
push!(k_c_fem_TM_2, 114689.36391759658)
push!(k_c_fem_TM_2, 114566.31611896376)
push!(k_c_fem_TM_2, 114526.44379572837)
push!(k_c_fem_TM_2, 114500.28550669279)


# Cargamos los resultados de orden 3
k_c_fem_TE_3 = Vector{Float64}()
push!(k_c_fem_TE_3, 18886.32072191986)
push!(k_c_fem_TE_3, 18886.317942531303)
push!(k_c_fem_TE_3, 18886.31785500316)
push!(k_c_fem_TE_3, 18886.31780042665)

k_c_fem_TM_3 = Vector{Float64}()
push!(k_c_fem_TM_3, 114500.83624038365)
push!(k_c_fem_TM_3, 114498.58906907425)
push!(k_c_fem_TM_3, 114498.37363664022)
push!(k_c_fem_TM_3, 114498.30314670774)



a_result = Numerical_integration.Mesh_object.Rect_waveguide_triangles.a
b_result = Numerical_integration.Mesh_object.Rect_waveguide_triangles.b

h = Vector{Float64}()
push!(h, a_result/5)
push!(h, a_result/8)
push!(h, a_result/10)
push!(h, a_result/20)

k_c_analytic_result = Numerical_integration.k_c_analytic

rel_error_result_1 = Vector{Float64}()
rel_error_result_2 = Vector{Float64}()
rel_error_result_3 = Vector{Float64}()

if Numerical_integration.FEM_TE == true
    for i in 1:length(h)
        error = abs((k_c_fem_TE_1[i]-k_c_analytic_result[1])/k_c_analytic_result[1])
        push!(rel_error_result_1, error)
    end

    for i in 1:length(h)
        error = abs((k_c_fem_TE_2[i]-k_c_analytic_result[1])/k_c_analytic_result[1])
        push!(rel_error_result_2, error)
    end

    for i in 1:length(h)
        error = abs((k_c_fem_TE_3[i]-k_c_analytic_result[1])/k_c_analytic_result[1])
        push!(rel_error_result_3, error)
    end
else
    for i in 1:length(h)
        error = abs((k_c_fem_TM_1[i]-k_c_analytic_result[3])/k_c_analytic_result[3])
        push!(rel_error_result_1, error)
    end
    
    for i in 1:length(h)
        error = abs((k_c_fem_TM_2[i]-k_c_analytic_result[3])/k_c_analytic_result[3])
        push!(rel_error_result_2, error)
    end

    for i in 1:length(h)
        error = abs((k_c_fem_TM_3[i]-k_c_analytic_result[3])/k_c_analytic_result[3])
        push!(rel_error_result_3, error)
    end
    
end


# Pintamos las graficas 

# Calcula el valor 1/h
x_values = 1 ./ h

# Grafica
plot(xaxis=:linear, yaxis=:log, legend=:bottomleft)
if Numerical_integration.FEM_TE == true
    plot!(x_values, rel_error_result_1, label="Orden 1", marker=:circle, linestyle=:solid)
    plot!(x_values, rel_error_result_2, label="Orden 2", marker=:circle, linestyle=:solid)
    plot!(x_values, rel_error_result_3, label="Orden 3", marker=:circle, linestyle=:solid)
    xlabel!("1/h")
    ylabel!("Relative Error")   
    title!("Modos TE - Error relativo vs Tamaño malla")
    ylims!(10e-15, 1.0)
    savefig("Results_images\\Mode_TE.png")
else
    plot!(x_values, rel_error_result_1, label="Orden 1", marker=:circle, linestyle=:solid)
    plot!(x_values, rel_error_result_2, label="Orden 2", marker=:circle, linestyle=:solid)
    plot!(x_values, rel_error_result_3, label="Orden 3", marker=:circle, linestyle=:solid)
    xlabel!("1/h")
    ylabel!("Relative Error")   
    title!("Modos TM - Error relativo vs Tamaño malla")
    ylims!(10e-15, 1.0)
    savefig("Results_images\\Mode_TM.png")
end

# Comprobar la pendiente para los diferentes nodos
if Numerical_integration.FEM_TE == true
    pen_1_TE = (log10(rel_error_result_1[2])-log10(rel_error_result_1[1]))/(log10(x_values[2])-log10(x_values[1]))
    pen_2_TE = (log10(rel_error_result_2[2])-log10(rel_error_result_2[1]))/(log10(x_values[2])-log10(x_values[1]))
    pen_3_TE = (log10(rel_error_result_3[2])-log10(rel_error_result_3[1]))/(log10(x_values[2])-log10(x_values[1]))
else
    pen_1_TM = (log10(rel_error_result_1[4])-log10(rel_error_result_1[3]))/(log10(x_values[4])-log10(x_values[3]))
    pen_2_TM = (log10(rel_error_result_2[4])-log10(rel_error_result_2[3]))/(log10(x_values[4])-log10(x_values[3]))
    pen_3_TM = (log10(rel_error_result_3[4])-log10(rel_error_result_3[3]))/(log10(x_values[4])-log10(x_values[3]))
end


