module NBodyIPs

using Reexport

# some generically useful code that
# could be used across different n-body basis function implementations
# TODO: move some codes from Invariants submodule to here
#       or maybe other parts of the package
include("misc.jl")

include("invariants.jl")
@reexport using NBodyIPs.Invariants

# describe basis functions in terms of symmetry invariants
include("polynomials.jl")
@reexport using NBodyIPs.Polynomials

# fitting from data (e.g., least squares)
include("fitting.jl")

# loading data
# TODO: move codes from the examples in here
# include("data.jl")

end # module
