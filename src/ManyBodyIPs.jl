module ManyBodyIPs

using Reexport

include("polynomials.jl")
# @reexport using Polynomials

include("assembly.jl")

include("fitting.jl")
# @reexport using Fitting

end # module
