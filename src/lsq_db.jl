
module DB

using JuLIP, ProgressMeter, Base.Threads, FileIO

using NBodyIPFitting: Dat, LsqDB
using NBodyIPFitting.Tools: tfor

import NBodyIPFitting.Lsq
const Lsq = NBodyIPFitting.Lsq

import Base: append!, push!

# export get_basis,
#        observations, get_lsq_system,
#        regularise, table,
#        config_types,
#        del_data!

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
   save(joinpath(dbdir, "data.jld2"), "data", Dat{Float64}[])
   save(joinpath(dbdir, "basis.jld2"), "basis", AbstractCalculator[])
   return LsqDB(dbdir)
end

function LsqDB(dbdir::AbstractString)
   if !isdir(dbdir)
      error("""`dirname` is not a directory. If you want to create a new LsqDB
               please use `initdb`.""")
   end
   data = load(joinpath(dbdir, "data.jld2"), "data")
   basis = load(joinpath(dbdir, "basis.jld2"), "basis")
   return LsqDB(data, basis, dbdir)
end

"""
`dbdir(db::LsqDB)` : return the absolute path to the database
"""
dbdir(db::LsqDB) = db.dirname

data(db::LsqDB) = db.data
basis(db::LsqDB) = db.basis

function push!(db::LsqDB, d::Dat;
               datidx = length(data(db)+1))
   datfname = joinpath(dbdir(db) * "dat_$(datidx).jld2")
   lsqdict = Lsq.evallsq(d, basis(db))
   save(datfname, "data", d, "lsq", lsqdict)
   # if length(db.data) >= datidx then we assume that
   # d has already been inserted into db.data
   if length(data(db)) < datidx
      @assert length(data(db)) == datidx-1
      push!(data(db), d)
   else
      @assert data(db)[datidx] == d
   end
   return db
end


function append!(db::LsqDB, ds::AbstractVector{TD}; verbose=true
                 ) where {TD <: Dat}
   # append the data
   len_data_old = length(data(db))
   append!(data(db), ds)
   # create the files in parallel
   tfor( n -> push!(db, ds[n]; datidx = len_data_old + n), 1:length(ds) )
end


function append!(db::LsqDB, b::AbstractVector{TB}) where {TB <: AbstractCalculator}
end



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
function _vec2arr(A_vec::AbstractVector{TA}) where {TA <: Array}
   B = vcat( [A[:] for A in A_vec]... )
   return reshape(B, tuple(size(A_vec[1])..., length(A_vec)))
end

_vec2arr(A_vec::AbstractVector{TA}) where {TA <: Real} = [a for a in A_vec]


   # sz = size(A_vec[1])
   # sz_arr = tuple(sz[1:end-1]..., length(A_vec))
   # A_arr = zeros(T, sz_arr)
   # colons = ntuple(_->:, N-1)
   # for i = 1:length(A_vec)
   #    A_arr[colons..., i] = A_vec[i]
   # end
   # return A_arr



end
