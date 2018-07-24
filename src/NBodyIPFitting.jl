
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

@reexport using JuLIP
@reexport using NBodyIPs

include("types.jl")

# loading data
include("data.jl")
@reexport using NBodyIPs.Data

include("lsq_db.jl")

# lsq and other fits
include("fitting.jl")

# assemble crude error tables and scatter plots
include("errors.jl")

#visualisation module
include("PIPplot.jl")
@reexport using NBodyIPs.PIPplot

end # module
