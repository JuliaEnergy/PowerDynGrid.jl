module DPSAGridTools

using Random
using LinearAlgebra
using SparseArrays

export build_nodal_admittance, line, trafo, Shunt, kron_reduction

begin # some useful helper functions
    # admittance from resistance and reactance
    yrx(r, x; pu=1) = (1 / (r + im * x)) / pu
    # admittance from conductance and susceptance
    ygb(g, b; pu=1) = (g + im * b) / pu
    # select only the offdiafonal part of a matrix
    # this should throw DimensionMismatch for non-square matrices
    maindiag(M) = spdiagm(0=>diag(M))
    offdiag(M) = M .- maindiag(M)
    # sparse inverse
    inv(M::AbstractSparseMatrix) = LinearAlgebra.inv(Array{Complex}(M))
end

# container of all possible grid elements
# TODO: refactor this to match PypSA conventions
abstract type GridElement end

# abstract types for grid elements (which might have different equivalent circuits) such that we can dispatch on them
abstract type OnePort <: GridElement end
abstract type TwoPort <: GridElement end
abstract type MultiPort <: GridElement end # three and more

struct Branch <: TwoPort
    nl::AbstractArray{Int} #node list
    circuit_matrix::AbstractArray # 2x2 line admittance matrix

    function Branch(k::Int, m::Int, y, y_shunt_km, y_shunt_mk, t_km, t_mk)
        """
        This branch representation uses the unified Π model.

        See also the Chapter 2 in
          Göran Andersson, _Power System Analysis_, Lecture 227-0526-00, ITET ETH Zurich, 2012

        Assumptions:
        * the line admittance is symmetric

        """
        Π = spzeros(Complex{Float64}, 2, 2)
        Π[1, 1] = abs2(t_km) * (y + y_shunt_km)
        Π[1, 2] = - conj(t_km) * t_mk * y
        Π[2, 1] = - conj(t_mk) * t_km * y
        Π[2, 2] = abs2(t_mk) * (y + y_shunt_mk)
        # Π[:, [k, m]] ./ pu # normalise to per unit admittance
        new([k, m], Π)
    end
end

# TODO: how can I dispatch the following as constructors of Branch?
trafo(lowV::Int, highV::Int, r, x, t_ratio) =
    Branch(lowV, highV, yrx(r, x), 0, 0, t_ratio, 1)
line(from::Int, to::Int, r, x, y_shunt) =
    Branch(from, to, yrx(r, x), y_shunt/2, y_shunt/2, 1, 1)

struct Shunt <: OnePort
    n::Int # node
    y_to_ground # admittance
    # TODO: check the sign, see p.17 of the ETH script
    Shunt(n, y) = new(n, -y)
    Shunt(n, r, x) = new(n, -yrx(r, x))
end

function add_element_admittance!(Y, el::TwoPort)
    Y[el.nl, el.nl] += el.circuit_matrix
    nothing
end

function add_element_admittance!(Y, el::OnePort)
    Y[el.n, el.n] += el.y_to_ground
    nothing
end

function build_nodal_admittance(number_of_nodes, element_list; verbose=false)
    Y = spzeros(Complex{Float64}, number_of_nodes, number_of_nodes)
    for el in element_list
        verbose ? println("add a ", typeof(el)) : nothing
        add_element_admittance!(Y, el)
    end
    Y
    #offdiag_sum = offdiag(Y) * ones(number_of_nodes)
    #maindiag(Y) + spdiagm(0 => offdiag_sum) - offdiag(Y) # return LY
end

function kron_reduction(Y, passive)
    Y[.~passive, .~passive] - Y[.~passive, passive] * inv(Y[passive, passive]) * Y[passive, .~passive]
    # TODO: does it make sense to return single values as a Number instead of Matrix?
    # size(Yred) == (1, 1) ? Yred[1, 1] : Yred
end

end # module
