
"""
# `NBodyIPFitting.jl`

Package for defining and fitting interatomic potentials based on the
N-Body expansion (ANOVA, HDMR, ...). The main exported type is
`NBodyIP` which is a `JuLIP` calculator.

See `?...` on how to
* `?NBodyIPFitting.Fitting` : fit an `NBodyIP`
* `?NBodyIPFitting.Data` : load data sets
* `?NBodyIPFitting.IO` : write and read an `NBodyIP`
"""
module NBodyIPFitting

using Reexport

include("tools.jl")

include("types.jl")

# loading data
include("data.jl")
@reexport using NBodyIPs.Data

include("lsq_db.jl")

# # lsq and other fits
# include("fitting.jl")
#
# # assemble crude error tables and scatter plots
# include("errors.jl")
#
# #visualisation module
# include("PIPplot.jl")
# @reexport using NBodyIPs.PIPplot

end # module
