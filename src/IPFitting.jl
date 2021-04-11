"""
# `IPFitting.jl`

Package for defining and fitting interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...). The main exported type is
`NBodyIP` which is a `JuLIP` calculator.

See `?...` on how to
* `?IPFitting.Fitting` : fit an `NBodyIP`
* `?IPFitting.Data` : load data sets
* `?IPFitting.IO` : write and read an `NBodyIP`
"""
__precompile__()
module IPFitting

using Reexport

include("tools.jl")

include("prototypes.jl")

include("obsiter.jl")

include("datatypes.jl")

# loading data
include("data.jl")
@reexport using IPFitting.Data

include("lsq_db.jl")
@reexport using IPFitting.DB

include("lsqerrors.jl")
@reexport using IPFitting.Errors

include("lsq.jl")
@reexport using IPFitting.Lsq

include("preconditioners.jl")
@reexport using IPFitting.Preconditioners

include("auxiliary.jl")

end # module
