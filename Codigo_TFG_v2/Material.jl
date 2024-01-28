module Material

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

end