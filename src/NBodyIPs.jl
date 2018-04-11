module NBodyIPs

using Reexport

# different standard dictionaries to use
include("dictionaries.jl")

# code for permutation invariant functions (usually polynomials of some sort)
include("polynomials.jl")

# code to combine the many-body terms from "polynomials.jl" into
# objects that can evaluate total energies and forces
include("assembly.jl")

# fitting from data (e.g., least squares)
include("fitting.jl")

# loading data
# include("data.jl")

end # module
