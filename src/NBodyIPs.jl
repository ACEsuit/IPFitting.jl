module NBodyIPs

export bodyorder

# prototypes
function bodyorder end

# miscallaneous
Base.ntuple(n::Integer, ::Val{M}) where {M} = ntuple(_->n, M)



# different standard dictionaries to use
include("dictionaries.jl")

# code for permutation invariant functions (usually polynomials of some sort)
include("polynomials.jl")

# describe basis functions in terms of symmetry invariants
include("invariants.jl")

# fitting from data (e.g., least squares)
include("fitting.jl")

# loading data
# include("data.jl")

end # module
