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
       INFO    KRON_________________ ..
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
    ___|___ __
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


using Base.Threads:          SpinLock
using StaticArrays:          SVector
using JuLIP:                 vecs, mat, AbstractCalculator, AbstractAtoms, Atoms
using NBodyIPFitting:        Dat, LsqDB, KronGroup, DataGroup, data, basis,
                             evaluate_lsq
using NBodyIPFitting.Data:   configtype, configname
using NBodyIPFitting.Tools:  tfor, decode
using NBodyIPs:              degree, bodyorder, basisname, combiscriptor

import NBodyIPFitting.FIO

import Base: flush, append!, union

export LsqDB, confignames

const KRONFILE = "_kron.h5"
const INFOFILE = "_info.json"

"""
`dbpath(db::LsqDB)` : return the absolute path to the database files, not
including the 'kron.h5' and 'info.json' endings.
"""
dbpath(db::LsqDB) = db.dbpath

infofile(dbpath::AbstractString) = dbpath * INFOFILE

kronfile(dbpath::AbstractString) = dbpath * KRONFILE

function load_info(dbpath::String)
   dbinfo = try
      FIO.load(infofile(dbpath))
   catch
      try
         FIO.load(dbpath * "_info.jld2")
      catch
         error("I tried to load both a json and a jld2 and neither worked?")
      end
   end
   basis = decode.(dbinfo["basis"])
   data_groups = Dict{String, DataGroup}()
   for (key, val) in dbinfo["data"]
      data_groups[key] = [Dat(v) for v in val]
   end
   return basis, data_groups
end

function save_info(dbpath::String, db)
   data = Dict()
   for (key, val) in db.data_groups
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
      # TODO: try to save it in a temp file
      return nothing
   end
   FIO.save(kronfile(dbpath), db.kron_groups)
end

save_kron(db) = save_kron(dbpath(db), db)

function load_kron(dbpath::String)
   kron_groups = FIO.load(kronfile(dbpath))
   return kron_groups
end

load_kron(db::LsqDB) = load_kron(dbpath(db))

function LsqDB(dbpath::AbstractString)
   basis, data_groups = load_info(dbpath)
   @assert isfile(kronfile(dbpath))
   kron_groups = load_kron(dbpath)
   return LsqDB(basis, data_groups, kron_groups, dbpath)
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
   # TODO: seems this fails to detect an existing database?
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
   config_types = unique(configtype.(data))
   # for key in config_types
   #    kron_groups[key] = KronGroup()
   # end
   db = LsqDB(basis, data_groups, kron_groups, dbpath)
   db_lock = SpinLock()
   lens = [length(d) for d in data]
   # parallel assembly of the LSQ matrix
   tfor( n -> safe_append!(db, db_lock, data[n]), 1:length(data);
         verbose=verbose, msg = "Assemble LSQ blocks",
         costs = lens)
   if dbpath != ""
      verbose && info("Writing db to disk...")
      try
         flush(db)
      catch
         warn("""something went wrong trying to save the db to disk, but the data
               should be ok; if it is crucial to keep it, try to save manually.""")
      end
      verbose && info("... done")
   else
      verbose && info("db is not written to disk since `dbpath` is empty.")
   end
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
      # add empty sub-groups for the datatypes (observationtypes) present in d
      for key in keys(d.D)
         D_ct[key] = Array{Float64}(length(d.D[key]), 0, length(db.basis))
      end
   else
      D_ct = db.kron_groups[ct]
      D_d = db.data_groups[ct]
   end
   # append d to INFO ...
   push!(D_d, d)
   # and append lsqrow to kron_groups[ct]
   # here, lsqrow is a dict with the same keys (=datatypes = observationtypes)
   # as kron_groups[ct], but containing 2D arrays
   for datatype in keys(d.D)
      db.kron_groups[ct][datatype] = _append(db.kron_groups[ct][datatype],
                                             lsqrow[datatype])
   end
   return nothing
end

function safe_append!(db::LsqDB, db_lock, d::Dat)
   # computing the lsq blocks ("rows") can be done in parallel,
   lsqrow = evallsq(d, basis(db))
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


function union(db1::LsqDB,db2::LsqDB; dbpath = (db1.dbpath * "_u"))
   info("Warning: the union implies that the data in the two databases are the same")
   configtypes = collect(keys(db1.data_groups))
   basis = cat(1, db1.basis, db2.basis)
   data_groups = db1.data_groups

   kron_groups1 = db1.kron_groups
   kron_groups2 = db2.kron_groups

   kron_groups = Dict{String, KronGroup}()
   for k in configtypes
      kron_groups[k] = KronGroup()
      @assert haskey(db2.data_groups, k)
      for ot in ["E", "F", "V"]
         if !haskey(data_groups[k][1].D, ot)
            continue
         end
         kron_groups[k][ot] = cat(3,kron_groups1[k][ot],kron_groups2[k][ot])
      end
   end
   db = LsqDB(basis, data_groups, kron_groups, dbpath)
   flush(db)
   return db
end

# function union(db1::LsqDB,db2::LsqDB, dbpath = (db1.dbpath * "_u"))
#    basis = [basis(db1), basis(db2)]
#    data_groups::Dict{String, DataGroup}
#       kron_groups::Dict{String, KronGroup}
#       dbpath::String
#    end
#
#    basis(db::LsqDB) = db.basis
#    basis(db::LsqDB, i::Integer) = db.basis[i]
#    data(db::LsqDB) = db.data
# end

# TODO: the following function suggests that the ordering of
#       <values, data, basis> was poorly chosen and it should instead be
#       <values, basis, data> => reconsider this!!!!

"""
if A = n x m x k and B = n x k, then this functions creates "appends"
B to A converting it into an n x (m+1) x k array.
 * n : the dimension of <data[i], basis[j]>
 * m : the number of Dat
 * k : the number of basis functions
"""
function _append(A::Array{T, 3}, B::Array{T, 2}) where {T}
   szA = size(A)
   if !(szA[1] == size(B,1) && szA[3] == size(B, 2))
      @show size(A), size(B)
      @assert (szA[1] == size(B,1) && szA[3] == size(B, 2))
   end
   Anew = Array{T}(szA[1], szA[2]+1, szA[3])
   Anew[:, 1:szA[2], :] .= A
   Anew[:, end, :] .= B
   return Anew
end


# ------------------- Evaluating LSQ Blocks -----------------


"""
Take a basis and split it into individual basis groups.
"""
function split_basis(basis)
   # get the types of the individual basis elements
   tps = combiscriptor.(basis)
   Iord = Vector{Int}[]
   Bord = Any[]
   for tp in unique(tps)
      # find which elements of basis have type `tp`
      I = find( [tp == t  for t in tps] )
      push!(Iord, I)
      push!(Bord, [b for b in basis[I]])
   end
   return Bord, Iord
end



# fill the LSQ system, i.e. evaluate basis at data points
function evallsq(d::Dat, B::AbstractVector{TB}) where {TB <: AbstractCalculator}
   B1 = [b for b in B]
   if !(isleaftype(eltype(B1)))
      return evallsq_split(d, B1)
   end
   # TB is a leaf-type so we can use "evaluate_many"
   return Dict( key => _evallsq(Val(Symbol(key)), B1, Atoms(d))
                for key in keys(d.D) )
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
   # vals will be a vector containing multiple evaluations
   vals = evaluate_lsq(vDT, B, at)
   # vectorise the first so we know the length of the data
   vec1 = vec(vDT, vals[1])
   # create a multi-D array to reshape these into
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
   A = Array{TA}(size(As[1],1), sum(size(AA, 2) for AA in As))
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
   sizes = Dict( key => [size(DD[key]) for DD in D_ord]
                 for key in keys(D_ord[1]) )
   return Dict( key => _cat_( [DD[key] for DD in D_ord], Iord )
                for key in keys(D_ord[1]) )
end


# ================ Convenience =============

confignames(db::LsqDB) = configname.( collect( keys( db.data_groups ) ) )

_nconfigs(db::LsqDB, configtype::AbstractString) =
      size( first(db.kron_groups[configtype])[2], 2 )

_nbasis(db::LsqDB, configtype::AbstractString) =
      size( first(db.kron_groups[configtype])[2], 3 )


function Base.info(db::LsqDB)
   # config names, how many
   configs_info = Dict{String, Int}()
   for (key, dg) in db.data_groups
      cn = configname(key)
      if haskey(configs_info, cn)
         configs_info[cn] += length(dg)
      else
         configs_info[cn] = length(dg)
      end
   end
   println("======================================================")
   println("       LsqDB Summary ")
   println("------------------------------------------------------")
   println(" Datagroup: configname  =>  number of configs ")
   for (i, (key, val)) in enumerate(configs_info)
      println("        $i : $key  =>  $val ")
   end
   # basis groups, how many, max degree
   B = db.basis
   Bord, Iord = split_basis(B)
   println("------------------------------------------------------")
   println(" Basis Group  |  type   | body-order |  degree | descriptor")
   for (i, B1) in enumerate(Bord)
      deg = maximum(degree.(B1))
      println("           $i  |  $(basisname(B1[1])) | $(bodyorder(B1[1])) | $deg |  " )
      # 
   end
   println("======================================================")
end

end
