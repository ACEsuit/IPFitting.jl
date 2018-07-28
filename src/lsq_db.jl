
module DB

using StaticArrays: SVector
using JuLIP: vecs, mat, AbstractCalculator

using NBodyIPFitting: Dat, LsqDB
using NBodyIPFitting.Tools: tfor, decode

import NBodyIPFitting.Lsq,
       NBodyIPFitting.FIO

import Base: append!, push!

# export get_basis,
#        observations, get_lsq_system,
#        regularise, table,
#        config_types,
#        del_data!

const DBFILETYPE = ".h5"
const DATAFILE = "data" * DBFILETYPE
const BASISFILE = "basis" * DBFILETYPE

# function _checkfname(fname)
#    if length(fname) >= 3
#       if fname[end-2:end] == ".h5"
#          return fname
#       end
#    end
#    return fname * ".h5"
# end

# function _save(fname, args...)
#    h5open(_checkfname(fname), "w") do file
#       for (key, val) in args
#          write(file, key, val)
#       end
#    end
# end
#
# function _load(fname, args...)
#    file = h5open(_checkfname(fname), "r")
#    ret = [file[key] for key in args]
#    close(file)
#    return ret...
# end


"""
`dbdir(db::LsqDB)` : return the absolute path to the database
"""
dbdir(db::LsqDB) = db.dirname

data(db::LsqDB) = db.data
data(db::LsqDB, i::Integer) = db.data[i]
basis(db::LsqDB) = db.basis
basis(db::LsqDB, i::Integer) = db.basis[i]

datafile(dbdir::AbstractString) = joinpath(dbdir, DATAFILE)
datafile(db::LsqDB) = datafile(dbdir(db))
load_data(dbdir::AbstractString) =
   Vector{Dat{Float64}}(Dat.(load(datafile(dbdir), "data")))
load_data(db::LsqDB) = load(dbdir(db))
save_data(dbdir::AbstractString, data) = _save(datafile(dbdir), "data" => Dict.(data))
save_data(db::LsqDB) = save_data(dbdir(db), data(db))

basisfile(dbdir::AbstractString) = joinpath(dbdir, BASISFILE)
basisfile(db::LsqDB) = basisfile(dbdir(db))
load_basis(dbdir::AbstractString) =
   Vector{AbstractCalculator}(decode.(load(basisfile(dbdir), "basis")))
load_basis(db::LsqDB) = load(dbdir(db))
save_basis(dbdir::AbstractString, basis) =
   _save(basisfile(dbdir), "basis" => Dict.(basis))
save_basis(db::LsqDB) = save_basis(dbdir(db), basis(db))

datfilename(db, datidx) = joinpath(dbdir(db), "dat_$(datidx)" * DBFILETYPE)

save_dat(db, datidx, lsq_dict) =
   _save(datfilename(db, datidx), "dat_" => data(db, datidx),
                                  collect(lsq_dict)...)

function load_dat(db, datidx)
   DD = _load(datfilename(db, datidx), "data", "lsq")
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
   save_data(dbdir, Dat{Float64}[])
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
   save(datfilename(db, datidx), "data", d, "lsq", lsqdict)
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

function append!(db::LsqDB, ds::AbstractVector{TD}; verbose=true) where {TD <: Dat}
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
   datfname = datfilename(db, datidx)
   lsqdict = load(datfname, "lsq")
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
      if key == "E"
         Darr[key] = val
      elseif key == "V"
         # each element of val is a R^6 vector representing a virial stress
         Darr[key] = mat(val)
      elseif key == "F"
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
      if key == "E"
         Dvec[key] = val
      elseif key == "V"
         Dvec[key] = reinterpret(SVector{6,Float64}, val, (size(val,2),))
      elseif key == "F"
         Dvec[key] = vecs.(_arr2vec(val))
      else
         error("unknown key `$key`")
      end
   end
   return Dvec
end


end
