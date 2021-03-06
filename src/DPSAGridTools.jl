module DPSAGridTools

using Random
using LinearAlgebra
using SparseArrays

export build_admittance, Line, Trafo, Shunt, kron_reduction

begin # some useful helper functions
    # admittance from resistance and reactance
    yrx(r, x; pu=1) = (1 / (r + im * x)) / pu
    # admittance from conductance and susceptance
    ygb(g, b; pu=1) = (g + im * b) / pu
    # select only the offdiafonal part of a matrix
    # this should throw DimensionMismatch for non-square matrices
    maindiag(M) = spdiagm(0=>diag(M))
    offdiag(M) = M .- maindiag(M)
    # to check if the resulting admittance is a Laplacian matrix
    islaplacian(M) = diag(M) + sum(offdiag(M), dims=2) == 0 ? true : false
    # sparse inverse
    inv(M::AbstractSparseMatrix) = LinearAlgebra.inv(Array{Complex}(M))
end

# container of all possible grid elements
# TODO: refactor this to match PypSA conventions
abstract type AbstractPort end

# abstract types for grid elements (which might have different equivalent circuits) such that we can dispatch on them
abstract type AbstractOnePort <: AbstractPort end
abstract type AbstractTwoPort <: AbstractPort end
abstract type AbstractMultiPort <: AbstractPort end # three and more

function pi_model(y, y_shunt_km, y_shunt_mk, t_km, t_mk)
    """
    Implementation of the unified Π model.

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
    Π
end

struct Line <: AbstractTwoPort
    nodes::AbstractArray #node list
    circuit_matrix::AbstractArray # 2x2 line admittance matrix

    function Line(from::Int, to::Int, r, x, y_shunt)
        """
        This line representation uses the Π model.

        """
        # Π[:, [k, m]] ./ pu # normalise to per unit admittance
        new([from, to], pi_model(yrx(r, x), y_shunt/2, y_shunt/2, 1, 1))
    end
end

struct Trafo <: AbstractTwoPort
    nodes::AbstractArray #node list
    circuit_matrix::AbstractArray # 2x2 line admittance matrix

    function Trafo(lowV::Int, highV::Int, r, x, t_ratio)
        """
        This transformer representation uses the Π model,
        assuming an ideal transformer in series with an admittance.
        The admittance is here taken to be on the high-voltage side.

        """
        # Π[:, [k, m]] ./ pu # normalise to per unit admittance
        new([lowV, highV], pi_model(yrx(r, x), 0, 0, t_ratio, 1))
    end
end

struct Shunt <: AbstractOnePort
    node # node
    y_to_ground # admittance
    # TODO: check the sign, see p.17 of the ETH script
    Shunt(n, y) = new(n, -y)
end

function add_element_admittance!(Y, el::AbstractTwoPort)
    Y[el.nodes, el.nodes] += el.circuit_matrix
    nothing
end

function add_element_admittance!(Y, el::AbstractOnePort)
    Y[el.node, el.node] += el.y_to_ground
    nothing
end

function build_admittance(number_of_nodes, element_list; verbose=false)
    Y = spzeros(Complex{Float64}, number_of_nodes, number_of_nodes)
    for el in element_list
        verbose ? println("add a ", typeof(el)) : nothing
        add_element_admittance!(Y, el)
    end
    Y
end

function kron_reduction(Y, passive)
    Y[.~passive, .~passive] - Y[.~passive, passive] * inv(Y[passive, passive]) * Y[passive, .~passive]
    # TODO: does it make sense to return single values as a Number instead of Matrix?
    # size(Yred) == (1, 1) ? Yred[1, 1] : Yred
end

end # module
