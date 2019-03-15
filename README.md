# NBodyIPFitting

<!-- [![Build Status](https://travis-ci.org/cortner/NBodyIPFitting.jl.svg?branch=master)](https://travis-ci.org/cortner/NBodyIPFitting.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/NBodyIPFitting.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/NBodyIPFitting.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/NBodyIPFitting.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/NBodyIPFitting.jl?branch=master) -->


## Basic Usage

### Step 1: Import data/observations

To import a database stored as an `xyz` file, use
```julia
data = NBodyIPFitting.Data.read_xyz(fname)
```
See `?read_xyz` for further options. This will return a `Vector{Dat}` where
each `Dat` is a container storing the atomistic configurion (`JuLIP.Atoms`),
the `configtype` as well as the "DFT observations".

### Step 2: Generate a basis

A basis is defined by
* choice of bond-length or bond-angle PIPs.
* space transform
* choice of cut-off
* body-order
* polynomial degree

For example, using bond-angle PIPs with Morse coordinates and a cosine cut-off
to model Si, we first define a descriptor
```julia
r0 = rnn(:Si)
rcut = 2.5 * r0
desc = BondAngleDesc("exp(- (r/$r0 - 1.0))", CosCut(rcut-1, rcut))
```
We can then generate basis functions using `nbpolys`, e.g.,
```julia
#            body-order, descriptor, degree
B4 = nbpolys(4,          desc,       8)
```
In practise, one would normally specify different cut-offs and space transforms
for different body-orders. Suppose these give descriptors `D2, D3, D4`, then
a 4-body basis can be constructed via
```julia
B = [ nbpolys(2, D2, 14); nbpolys(3, D3, 11); nbpolys(4, D4, 8) ]
```

For more details and more complex basis sets, see below.


### Step 3: Precompute a Lsq system

Once the dft dataset and the basis functions have been specified, the
least-squares system matrix can be assembled. This can be very time-consuming
for high body-orders, large basis sets and large data sets. Therefore this
matrix is stored in a block format that allows us to later re-use it in a variety
of different ways. This is done via
```julia
db = LsqDB(fname, configs, basis)
```
* The db is stored in two files: `fname_info.jld2` and `fname_kron.h5`. In
particular, `fname` is the path + name of the db, but without ending. E.g,
`"~/scratch/nbodyips/W5Bdeg12env"`.
* `configs` is a `Vector{Dat}`
* `basis` is a `Vector{<: AbstractCalculator}`
* The command `db = LsqDB(fname, configs, basis)` evaluates the basis functions,
e.g.,  `energy(basis[i], configs[j].at)` for all `i, j`, and stores these values
which make up the lsq matrix.

To reload a pre-computed lsq system, use `LsqDB(fname)`. To compute a lsq
system without storing it on disk, use `LsqDB("", configs, basis)`, i.e.,
pass an empty string as the filename.

### Step 4: Lsq fit, Analyse the fitting errors

The main function to call is
`lsqfit(db; kwargs...) -> IP, fitinfo`
The system is solved via (variants of) the QR factorisation. See `?lsqfit`
for details.

### Step 5: Usage

The output `IP` of `lsqfit` is a `JuLIP.AbstractCalculator` which supports
`energy, forces, virial, site_energies`. (todo: write more here, in
particular mention `fast`)


## More on Descriptors and Basis sets

### Bond-length descriptor

### Bond-angle descriptor

### Environment-dependent N-Body potentials


## Regularisation


## Hooks
