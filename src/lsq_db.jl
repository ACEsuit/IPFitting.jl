"""
`module DB`

This module implements a very basic "database" for storing a precomputed
LSQ system. This is useful for fitting NBodyIPs since it allows one to
precompute the <basis, data> inner products which are very expensive,
and then quickly construct LSQ systems from them using e.g., many variants
of weighting and regularisation.

A db called (e.g.) 'lsqdb' consists of two files:
 * `lsqdb_info.json` : this stores a list of basis functions and a list of
"data", i.e. configurations
 * `lsqdb_kron.h5` : this stores all the inner products <basis, data>

## The INFO file

### Basis Functions

The basis is (for the time being) simply represented as a list (Vector) of
JuLIP.AbstractCalculator. These are stored as `INFO["basis"]`.

### Data Groups

Each "datum" is an atomistic configuration, which must have a "configtype".
Two data with the same configtype are required to contain the same number of
atoms (informally, they should also represent "similar" kinds of configurations).
`data = INFO["data"]` is a Dictionary where the keys are the configtypes.
E.g. suppose there is a configtype "vacancy", then `d = data["vacancy"]` is
a vector of `Dat` instances (stored as Dictionaries); see `?Dat` for more
information.


## The KRON file

The KRON file contains several dictionaries (groups in HDF5 terminology)
where each group represents a configtype. E.g., take a configtype "vacancy",
then "D = KRON["vacancy"] is a dictionary where each key represents a type of
data, e.g., energy, forces, virial. E.g., let the key "F" represent forces,
then `d["F"]` is a 3-dimensional array where `d["F"][i, j, :]` contains the
forces obtained from <basis[i], data[j]> into vectorised format.

### (De-)Vectorising Data

Forces in `JuLIP` are represented as `Vector{SVector{...}}`, which is the
same memory layout as a 3 x N matrix (`JuLIP.mat`), which can can then be
vectorised (`[:]`), and this vectoriation is readily undone again.
Analogously, any data must be stored in such a vectorised format. This is
achieved (e.g. for forces) via
* `Base.vec(::Val{:F}, F) -> Vector`
* `NBodyIPFitting.devec(::Val{:F}, Fvec) -> Vector{SVector{...}}`

## Usage

"""
module DB

using StaticArrays: SVector
using JuLIP: vecs, mat, AbstractCalculator
using FileIO

using NBodyIPFitting: Dat, LsqDB
using NBodyIPFitting.Tools: tfor, decode

import NBodyIPFitting.Lsq,
       NBodyIPFitting.FIO

import Base: append!, push!


const KRONFILE = "_kron.h5"
const INFOFILE = "_info.json"


"""
`dbdir(db::LsqDB)` : return the absolute path to the database
`dbpath(db::LsqDB)` : return the absolute path to the database
"""
dbdir(db::LsqDB) = db.dirname
dbpath(db::LsqDB) = db.dirname

data(db::LsqDB) = db.data
data(db::LsqDB, i::Integer) = db.data[i]
basis(db::LsqDB) = db.basis
basis(db::LsqDB, i::Integer) = db.basis[i]

datafile(dbdir::AbstractString) = joinpath(dbdir, DATAFILE)
datafile(db::LsqDB) = datafile(dbdir(db))
load_data(dbdir::AbstractString) =
   Vector{Dat}(Dat.(load(datafile(dbdir), "data")))
load_data(db::LsqDB) = load(dbdir(db))
save_data(dbdir::AbstractString, data) = save(datafile(dbdir), "data", Dict.(data))
save_data(db::LsqDB) = save_data(dbdir(db), data(db))

basisfile(dbdir::AbstractString) = joinpath(dbdir, BASISFILE)
basisfile(db::LsqDB) = basisfile(dbdir(db))
load_basis(dbdir::AbstractString) =
   Vector{AbstractCalculator}(decode.(load(basisfile(dbdir), "basis")))
load_basis(db::LsqDB) = load(dbdir(db))
save_basis(dbdir::AbstractString, basis) = save(basisfile(dbdir), "basis", Dict.(basis))
save_basis(db::LsqDB) = save_basis(dbdir(db), basis(db))

datfilename(db, datidx) = joinpath(dbdir(db), "dat_$(datidx)" * DBFILETYPE)

save_dat(db, datidx, lsq_dict) =
   FIO.save(datfilename(db, datidx), "dat" => Dict(data(db, datidx)),
                                     collect(_vec2arr(lsq_dict))...)

function load_dat(db, datidx)
   DD = _load(datfilename(db, datidx), "dat", "lsq")
   @assert dat == data(db, datidx)
   return _arr2vec(lsq)
end

"""
`function initdb(basedir, dbname)`

Initialise an empty lsq database.

* `basedir`: base directory where the database will be stored
* `dbname`: name of the database (this will be the name of a directory where
all files will be stored
"""
function initdb(basedir, dbname)
   @assert !('/' in dbname)
   @assert isdir(basedir)
   dbdir = joinpath(basedir, dbname)
   @assert !isdir(dbdir)
   mkdir(dbdir)
   save_data(dbdir, Dat[])
   save_basis(dbdir, AbstractCalculator[])
   return LsqDB(dbdir)
end

function LsqDB(dbdir::AbstractString)
   if !isdir(dbdir)
      error("""`dbdir` is not a directory. If you want to create a new `LsqDB`
               please use `initdb`.""")
   end
   return LsqDB(load_data(dbdir), load_basis(dbdir), dbdir)
end

# ------------- Append New Data to the DB -----------------

function push!(db::LsqDB, d::Dat; datidx = length(data(db)+1))
   # TODO: check whether d already exists in the database
   lsqdict = Lsq.evallsq(d, basis(db))
   save_dat(db, datidx, lsqdict)
   # if length(db.data) >= datidx then we assume that d has already been
   # inserted into db.data, otherwise push it
   if length(data(db)) < datidx
      @assert length(data(db)) == datidx-1
      push!(data(db), d)
      save_data(db)
   else
      @assert data(db)[datidx] == d
   end
   return db
end

function append!(db::LsqDB, ds::AbstractVector{Dat}; verbose=true)
   # append the data
   len_data_old = length(data(db))
   append!(data(db), ds)
   save_data(db)
   # create the dat_x_basis files in parallel
   tfor( n -> push!(db, ds[n]; datidx = len_data_old + n), 1:length(ds);
         verbose=verbose, msg = "Assemble LSQ blocks")
   return db
end

# ------------- Append New Basis Functions to the DB -----------------

push!(db::LsqDB, b::AbstractCalculator; kwargs...) = append!(db, [b]; kwargs...)

function append!(db::LsqDB, bs::AbstractVector{TB};
                 verbose=true) where {TB <: AbstractCalculator}
   # TODO : check whether bs already exist in the database?
   append!(basis(db), bs)
   save_basis(db)
   # loop through all existing datafiles and append the new basis functions
   # to each of them
   tfor( n -> _append_basis_to_dat!(db, bs, n), 1:length(data(db));
         verbose=verbose, msg="Append New LSQ blocks" )
   return db
end

function _append_basis_to_dat!(db, bs, datidx)
   lsqdict = load_dat(lsq, datidx)
   lsqnew = Lsq.evallsq(data(db, datidx), bs)
   for key in keys(lsqdict)
      append!(lsqdict[key], lsqnew[key])
   end
   save_dat(db, datidx, lsqdict)
   return db
end

# --------- Auxiliaries -------------


"""
`function _vec2arr(As)`

concatenate several arrays along the last dimension
"""
function _vec2arr(A_vec::AbstractVector{TA}) where {TA <: Array}
   B = vcat( [A[:] for A in A_vec]... )
   return reshape(B, tuple(size(A_vec[1])..., length(A_vec)))
end

_vec2arr(A_vec::AbstractVector{TA}) where {TA <: Real} = [a for a in A_vec]

"""
`function _arr2vec(A_arr)`

split a long multi-dimensional array into a vector of lower-dimensional
arrays along the last dimension.
"""
function _arr2vec(A_arr::AbstractArray{TA}) where {TA <: Real}
   colons = ntuple(_->:, length(size(A_arr))-1)
   return [ A_arr[colons..., n] for n = 1:size(A_arr)[end] ]
end

_arr2vec(A_vec::AbstractVector{TA}) where {TA <: Real} = [a for a in A_vec]


"""
`function _vec2arr(Dvec::Dict{String})`

convert LSQ matrix entries stored in atomistic/JuLIP format into a
big multi-dimensional arrays for efficient storage
"""
function _vec2arr(Dvec::Dict{String})
   Darr = Dict{String, Array}()
   for (key, val) in Dvec
      if key == ENERGY
         Darr[key] = val
      elseif key == VIRIAL
         # each element of val is a R^6 vector representing a virial stress
         Darr[key] = mat(val)
      elseif key == FORCES
         # mat.(val) : convert each collection of forces into a matrix
         # _vec2arr( : convert a Vector of matrices into a multi-dimensional array
         Darr[key] = _vec2arr( mat.(val) )
      else
         error("unknown key `$key`")
      end
   end
   return Darr
end

"""
function _arr2vec(Darr::Dict{String})`

convert LSQ matrix entries stored as multi-dimensional arrays back into
atomistic/JuLIP format
"""
function _arr2vec(Darr::Dict{String,Array})
   Dvec = Dict{String, Vector}()
   for (key, val) in Darr
      if key == ENERGY
         Dvec[key] = val
      elseif key == VIRIAL
         Dvec[key] = reinterpret(SVector{6,Float64}, val, (size(val,2),))
      elseif key == FORCES
         Dvec[key] = vecs.(_arr2vec(val))
      else
         error("unknown key `$key`")
      end
   end
   return Dvec
end


end
