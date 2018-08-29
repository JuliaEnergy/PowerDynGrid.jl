module PowerDynGridTools

using Random
using LinearAlgebra
using SparseArrays

begin # some useful helper functions
    # admittance from resistance and reactance
    yrx(r, x; pu=1) = (1 / (r + im * x)) / pu
    # admittance from conductance and susceptance
    ygb(g, b; pu=1) = (g + im * b) / pu
end

struct PiModel
    node_range
    circuit_matrix
    incidence_matrix
end

function pi_model(k, m, y, y_shunt_km, y_shunt_mk, t_km, t_mk)
    # See also the Chapter 2 in
    #   Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012
    #
    # Assumptions:
    # * the line admittance is symmetric
    B = spzeros(Complex, 2, 2)
    B[1, 1] = abs2(t_km) * (y + y_shunt_km)
    B[1, 2] = - conj(t_km) * t_mk * y
    B[2, 1] = - conj(t_mk) * t_km * y
    B[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
    K = spdiagm(0=>[1, 1])
    PiModel(k:m, B, K)
end
# B[:, [k, m]] ./ pu # normalise to per unit admittance

function trafo_twoway(lowV, highV, r, x, t_ratio)
    # admittance is located at the high-voltage end
    pi_model(lowV, highV, yrx(r, x), 0, 0, t_ratio, 1)
end

function line(from, to, r, x, y_shunt)
    # admittance is located at the high-voltage end
    pi_model(from, to, yrx(r, x), y_shunt, y_shunt, 1, 1)
end


# export pi_model

end # module
