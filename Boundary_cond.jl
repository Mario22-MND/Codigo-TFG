module Boundary_cond

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

end