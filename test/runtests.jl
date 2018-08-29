using PowerDynGridTools

using Test
using SparseArrays

nodes = 3
edges = 2

B = spzeros(Complex, 2edges, nodes)
K = spzeros(Int, 2edges, nodes)

for i in 1:nodes-1
    branch = PowerDynGridTools.line(i, i+1, i*0.0062, 0.036, 0.0105im)
    B[i:i+1, branch.node_range] = branch.circuit_matrix
    K[i:i+1, branch.node_range] = branch.incidence_matrix
end

@show Y = K' * B
