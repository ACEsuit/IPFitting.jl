module ManyBodyIPs

# package code goes here

include("polynomials.jl")
@reexport using Polynomials

include("assembly.jl")

include("fitting.jl")
@reexport using Fitting

end # module
