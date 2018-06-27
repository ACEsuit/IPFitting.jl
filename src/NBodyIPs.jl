
"""
# `NBodyIPs.jl`

Package for defining and fitting interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...). The main exported type is
`NBodyIP` which is a `JuLIP` calculator.

See `?...` on how to
* `?NBodyIPs.Polys` : specify a polynomial basis set
* `?NBodyIPs.Fitting` : fit an `NBodyIP`
* `?NBodyIPs.Data` : load data sets
* `?NBodyIPs.IO` : write and read an `NBodyIP`
"""
module NBodyIPs

using Reexport

@reexport using StaticArrays
@reexport using FileIO
@reexport using JuLIP


# two auxiliary functions to make for easier assembly of the code
# TODO: move these somewhere else
push_str!(ex::Vector{Expr}, s::String) = push!(ex, parse(s))
append_str!(ex::Vector{Expr}, s::Vector{String}) = append!(ex, parse.(s))


include("fastpolys.jl")

include("invariants.jl")

include("common.jl")


# some generically useful code that
# could be used across different n-body basis function implementations
# TODO: move some codes from Invariants submodule to here
#       or maybe other parts of the package
include("misc.jl")


# describe basis functions in terms of symmetry invariants
include("polynomials.jl")
@reexport using NBodyIPs.Polys

# loading data
include("data.jl")
@reexport using NBodyIPs.Data

# fitting from data (e.g., least squares)
# TODO: make this a sub-module and re-export
include("fitting.jl")

include("errors.jl")

# # IP i/o
# include("io.jl")

#visualisation module
include("PIPplot.jl")
@reexport using NBodyIPs.PIPplot

end # module
