"""
`module DB`

## Overview

This module implements a very basic "database" for storing a precomputed
LSQ system. This is useful for fitting NBodyIPs since it allows one to
precompute the <basis, data> inner products which are very expensive,
and then quickly construct LSQ systems from them using e.g., many variants
of weighting and regularisation.

A db called (e.g.) 'lsqdb' consists of two files:
 * `lsqdb_info.xxxx` : (`xxxx` will be either `json` or `jld2` - to be decided)
this stores a list of "basis" functions and a list of "data", i.e. configurations
 * `lsqdb_kron.h5` : this stores all the inner products <basis, data>

            DB
           /  \
       INFO    KRON_________________ ...
     /     \           |    |    |
  BASIS    DATA       CT1  CT2  CT3
       _____|_____ ...
      |    |    |
     CT1  CT2  CT3

 * BASIS : Vector{AbstractCalculator}
 * DATA : Vector{Dat}
 * CTj < INFO: config types; a `Vector{Dat}`;
 * CTj < KRON: a dictionary with the following structure

       CT1
    ___|___ __...
   |   |   |
   E   F   V

where E, F, V are Array{Float64, 3}, these ids can in principle by anything
but need to be "registered in `datatypes.jl`

## More details and remarks

### INFO

* BASIS: The basis is (for the time being) simply represented as a list (Vector) of
JuLIP.AbstractCalculator. These are stored as `INFO["basis"]`.

* DATA: Each "datum" is an atomistic configuration, which must have a "configtype".
Two data with the same configtype are required to contain the same number of
atoms (informally, they should also represent "similar" kinds of configurations).
`data = INFO["data"]` is a Dictionary where the keys are the configtypes.
E.g. suppose there is a configtype "vacancy", then `d = data["vacancy"]` is
a vector of `Dat` instances (stored as Dictionaries); see `?Dat` for more
information.

### KRON

The KRON file contains several dictionaries (groups in HDF5 terminology)
where each group represents a configtype. E.g., take a configtype "vacancy",
then "D = KRON["vacancy"] is a dictionary where each key represents a type of
data, e.g., energy, forces, virial. E.g., let the key "F" represent forces,
then `d["F"]` is a 3-dimensional array where `d["F"][:, i, j]` contains the
forces obtained from <data[i], basis[j]> into vectorised format.

### (De-)Vectorising Data

Forces in `JuLIP` are represented as `Vector{SVector{...}}`, which is the
same memory layout as a 3 x N matrix (`JuLIP.mat`), which can can then be
vectorised (`[:]`), and this vectoriation is readily undone again.
Analogously, any data must be stored in such a vectorised format. This is
achieved (e.g. for forces) via
* `Base.vec(::Val{:F}, F) -> Vector`
* `NBodyIPFitting.devec(::Val{:F}, Fvec) -> Vector{SVector{...}}`

## Usage and Examples

"""
module DB


using Threads:               SpinLock
using StaticArrays:          SVector
using JuLIP:                 vecs, mat, AbstractCalculator
using NBodyIPFitting:        Dat, LsqDB, KronGroup, DataGroup
using NBodyIPFitting.Tools:  tfor, decode

import NBodyIPFitting.Lsq,
       NBodyIPFitting.FIO

import Base: flush, append!

const KRONFILE = "_kron.h5"
const INFOFILE = "_info.jld2"

"""
`dbpath(db::LsqDB)` : return the absolute path to the database files, not
including the '_kron.h5' and '_info.json' endings.
"""
dbpath(db::LsqDB) = db.dbpath

basis(db::LsqDB) = db.basis
basis(db::LsqDB, i::Integer) = db.basis[i]

data(db::LsqDB) = db.data
# data(db::LsqDB, i::Integer) = db.data[i]

infofile(dbpath::AbstractString) = dbpath * INFOFILE

kronfile(dbpath::AbstractString) = dbpath * KRONFILE

function load_info(dbpath)
   info = FIO.load(infofile(dbpath))
   basis = decode.(info["basis"])
   data_groups = Dict{String, Dat}()
   for (key, val) in info["data"]
      data_groups[key] = Dat.(val)
   end
   return basis, data_groups
end

function save_info(dbpath, db)
   data = Dict()
   for (key, val) in db.data
      data[key] = Dict.(val)
   end
   FIO.save(infofile(dbpath),
            Dict("basis" => Dict.(db.basis), "data" => data))
   return nothing
end

save_info(db) = save_info(dbpath(db), db)

function save_kron(dbpath, db)
   if isfile(kronfile(dbpath))
      warn("""trying to save `kron`, but kronfile already exists; aborting""")
   end
   FIO.save(kronfile(dbpath), db.kron_groups)
end

save_kron(db) = save_kron(dbpath(db), db)

function LsqDB(dbpath::AbstractString)
   basis, data_groups = load_info(dbpath)
   @assert isfile(kronfile(dbpath))
   return LsqDB(basis, data_groups, Dict{String, KronGroup}(), dbpath)
end


"""
for the time being, this just checks that the db director and db name are
admissible. In the future, this should initialise the DB files.
"""
function initdb(basedir, dbname)
   @assert !('/' in dbname)
   @assert isdir(basedir)
   dbpath = joinpath(basedir, dbname)
   initdb(dbpath)
   return nothing
end

function initdb(dbpath)
   @assert !isfile(infofile(dbpath))
   @assert !isfile(kronfile(dbpath))
   # check that we can actually create and delete this file
   FIO.save(infofile(dbpath), Dict("basis" => [], "data" => []))
   rm(infofile(dbpath))
   FIO.save(kronfile(dbpath), Dict("a" => rand(10), "b" => rand()))
   rm(kronfile(dbpath))
   return nothing
end



function LsqDB(dbpath::AbstractString,
               basis::AbstractVector{<: AbstractCalculator},
               data::AbstractVector{Dat};
               verbose=true)
   data_groups = Dict{String, DataGroup}()
   kron_groups = Dict{String, KronGroup}()
   config_types = unique(config_type.(data))
   for key in config_types
      kron_groups[key] = KronGroup()
   end
   db = LsqDB(basis, data_groups, kron_groups, dbpath)
   db_lock = SpinLock()
   # parallel assembly of the LSQ matrix
   tfor( n -> safe_append!(db, db_lock, data[n]), 1:length(data);
         verbose=verbose, msg = "Assemble LSQ blocks" )
   verbose && info("Writing db to disk...")
   try
      flush(db)
   catch
      warn("""something went wrong trying to save the db to disk, but the data
            should be ok; if it is crucial to keep it, try to save manually.""")
   end
   verbose && info("... done")
   return db
end


function _append_inner!(db::LsqDB, d::Dat, lsqrow::Dict)
   ct = configtype(d)
   # if this config_type does not yet exist in db.kron_groups, then
   # initialise this group.
   if !haskey(db.kron_groups, ct)
      @assert !haskey(db.data_groups, ct)
      D_ct = db.kron_groups[ct] = KronGroup()
      D_d = db.data_groups[ct] = DataGroup()
      # add empty sub-groups for the datatypes present in d
      for key in keys(d.D)
         D_ct[key] = Array{Float64, 3}()
      end
   end
   # append d to INFO ...
   push!(D_d[ct], d)
   # and append lsqrow to kron_groups[ct]
   # here, lsqrow is a dict with the same keys as kron_groups[ct], but
   # containing 2D arrays
   for datatype in keys(d.D)
      db.kron_groups[ct][datatype] = _append!(db.kron_groups[ct][datatype],
                                              lsqrow[datatype])
   end
   return nothing
end

function safe_append!(db::LsqDB, db_lock, d::Dat)
   # computing the lsq blocks ("rows") can be done in parallel,
   lsqrow = Lsq.evallsq(d, basis(db))
   # but writing them into the DB must be done in a threadsafe way
   lock(db_lock)
   _append_inner!(db, d, lsqrow)
   unlock(db_lock)
   return nothing
end

function flush(db::LsqDB)
   save_info(db)
   save_kron(db)
   return nothing
end

"""
if A = n x m x k and B = n x m, then this functions creates "appends"
B to A converting it into an n x m x (k+1) array.
"""
function _append!(A::Array{T, 3}, B::Array{T, 2}) where {T}
   szA = size(A)
   @assert szA[1:2] == size(B)
   Avec = vec(A)
   append!(Avec, B[:])
   return reshape(Avec, szA[1], szA[2], szA[3]+1)
end


# ------------------- Evaluating LSQ Blocks -----------------

"""
Take a basis and split it into individual basis groups.
"""
function split_basis(basis)
   # get the types of the individual basis elements
   tps = typeof.(basis)
   Iord = Vector{Int}[]
   Bord = Any[]
   for tp in unique(tps)
      # find which elements of basis have type `tp`
      I = find( tp .== tps )
      push!(Iord, I)
      push!(Bord, [b for b in basis[I]])
   end
   return Bord, Iord
end


# fill the LSQ system, i.e. evaluate basis at data points
function evallsq(d::Dat, B::AbstractVector{TB}
                 ) where {TB <: AbstractCalculator}
   if !(isleaftype(TB))
      return evallsq_split(d, B)
   end
   # TB is a leaf-type so we can use "evaluate_many"
   return Dict( key => _evallsq(Val(Symbol(key)), B, Atoms(d)) )
end

"""
evaluate one specific kind of datum, such as energy, forces, etc
for one Dat across a large section of the basis; implicitly this will
normally be done via `evaluate_many` or similar. The line
```
vals = evaluate_lsq(vDT, B, at)
```
should likely be the bottleneck of `evallsq`
"""
function _evallsq(vDT::Val,
                  B::AbstractVector{<:AbstractCalculator},
                  at::AbstractAtoms)
   vals = evaluate_lsq(vDT, B, at)
   vec1 = vec(vDT, vals[1])
   A = Array{Float64}(length(vec1), length(B))
   A[:, 1] = vec1
   for n = 2:length(B)
      A[:, n] = vec(vDT, vals[n])
   end
   return A
end

"""
concatenate several arrays Ai of shape ni x m  along the last dimension into
an array ∑ni x m.
"""
function _cat_(As, Iord)
   TA = eltype(As[1])
   A = Vector{TA}(size(As[1],2),
                  sum(size(AA, 1) for AA in As))
   for i = 1:length(As)
      A[:, Iord[i]] = As[i]
   end
   return A
end

"""
split the Basis `B` into subsets of identical types and evaluate
those independently (fast).
"""
function evallsq_split(d, basis)
   # TB is not a leaf-type so we should split the basis to be able to
   # evaluate_many & co
   Bord, Iord = split_basis(basis)
   D_ord = [ evallsq(d, BB)  for BB in Bord ]
   # each D_ord[i] is a Dict storing the LSQ system components.
   # we assume that all D_ord[i] contain the same keys.
   return Dict( key => _cat_( [DD[key] for DD in D_ord], Iord ),
                for key in keys(D_ord[1]) )
end


end
