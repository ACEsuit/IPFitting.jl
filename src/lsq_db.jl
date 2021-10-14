"""
`module DB`
## Overview
This module implements a "database" for storing a precomputed LSQ system. This
is useful for fitting e.g. NBodyIPs since it allows one to precompute the <basis,
data> inner products which are very expensive, and then quickly construct LSQ
systems from them using e.g., many variants of weighting and regularisation.
A db called (e.g.) 'lsqdb' consists of two files:
 * `lsqdb_info.json` : this stores a list of basis functions and a list
 of configurations with attached observations (normally DFT data)
 * `lsqdb_kron.h5` : this stores all the inner products <basis, data>
 as a single large matrix. The observations and basis functions in
 `lsqdb_info` store the corresponding column and row indices.
## More details and remarks
### INFO
* BASIS: The basis is (for the time being) simply represented as a list (Vector) of
JuLIP.AbstractCalculator. These are stored as `INFO["basis"]`.
* DATA: Each "datum" is an atomistic configuration, which must have a
"configtype", stored as a `Dat`.
Informally two pieces of data with the same configtype should represent
"similar" kinds of configurations. `data = INFO["data"]` is a Dictionary
where the keys are the configtypes.
### (De-)Vectorising Data
Forces in `JuLIP` are represented as `Vector{SVector{...}}`, which is the
same memory layout as a 3 x N matrix (`JuLIP.mat`), which can then be
vectorised (`[:]`), and this vectoriation is readily undone again.
Analogously, any data must be stored in such a vectorised format. This is
achieved (e.g. for forces) via
* `vec_obs(::Val{:F}, F) -> Vector`
* `devec_obs(::Val{:F}, Fvec) -> Vector{SVector{...}}`
or equivalently
* `vec_obs("F", F) -> Vector`
* `devec_obs("F", Fvec) -> Vector{SVector{...}}`
(the `Val` versions are used for performance optimisation)
"""
module DB

using Base.Threads:          SpinLock, nthreads
using StaticArrays:          SVector
using JuLIP:                 AbstractCalculator, AbstractAtoms, Atoms, energy, forces,
                             save_dict, load_dict, read_dict, write_dict
using JuLIP.MLIPs:           IPBasis
using IPFitting:        Dat, LsqDB, basis, eval_obs, observations,
                             observation, vec_obs, devec_obs,
                             tfor_observations
using IPFitting.Data:   configtype
using HDF5:                  h5open, read

using DistributedArrays
using Distributed

import Base: flush, append!, union

export LsqDB, LsqDB_dist, info, configtypes

const VERSION = 2
const KRONFILE = "_kron.h5"
const INFOFILE = "_info.json"

"""
`dbpath(db::LsqDB)` : return the absolute path to the database files, not
including the 'kron.h5' and 'info.json' endings.
"""
dbpath(db::LsqDB) = db.dbpath

infofile(dbpath::AbstractString) = dbpath * INFOFILE

kronfile(dbpath::AbstractString) = dbpath * KRONFILE

# ------------ Save and load the info file

"if the file exists, append a random string to avoid overwriting"
function _backupfile(fname)
   if isfile(fname)
      fnew = fname * "." * String(rand('a':'z', 5))
      @warn("The file $fname already exists. It will be renamed to $fnew to avoid overwriting.")
      run(`mv $fname $fnew`)
   end
   return nothing
end

function load_info(dbpath::String)
   dbinfo = load_dict(infofile(dbpath))
   if !haskey(dbinfo, "version")
      @error("This LsqDB was produced with an older version of IPFitting.")
   end
   version = dbinfo["version"]
   if version < VERSION
      @error("This LsqDB was produced with an older version of IPFitting.")
   end
   basis = read_dict(dbinfo["basis"])
   configs = Dat.(dbinfo["configs"])   # here we already know the type
   return basis, configs
end

function save_info(dbpath::String, db)
   _backupfile(infofile(dbpath))
   save_dict(infofile(dbpath),
             Dict("version" => VERSION,
                  "basis" => write_dict(db.basis),
                  "configs" => write_dict.(db.configs))
            )
   return nothing
end

save_info(db) = save_info(dbpath(db), db)

# ----------- Save and load the KRON file

"save a single matrix to HDF5"
_savemath5(A, fname) =
   h5open(fname, "w") do fid
      fid["A"] = A
      nothing
   end

"load a single matrix from HDF5"
_loadmath5(fname) =
   h5open(fname, "r") do fid
      read(fid["A"])
   end

function save_kron(dbpath, db)
   _backupfile(kronfile(dbpath))
   _savemath5(db.Ψ, kronfile(dbpath))
end

save_kron(db) = save_kron(dbpath(db), db)

load_kron(dbpath::String; mmap=false) = _loadmath5(kronfile(dbpath))

load_kron(db::LsqDB; mmap=false) = load_kron(dbpath(db); mmap=mmap)

function LsqDB(dbpath::AbstractString; mmap=true)
   basis, configs = load_info(dbpath)
   Ψ = load_kron(dbpath; mmap=mmap)
   return LsqDB(basis, configs, Ψ, dbpath)
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
   save_dict(infofile(dbpath), Dict("version" => VERSION, "basis" => [], "data" => []))
   rm(infofile(dbpath))
   _savemath5(rand(5,5), kronfile(dbpath))
   rm(kronfile(dbpath))
   return nothing
end

function flush(db::LsqDB)
   save_info(db)
   save_kron(db)
   return nothing
end

function LsqDB(dbpath::AbstractString,
               basis::IPBasis,
               configs::AbstractVector{Dat};
               verbose=true,
               maxnthreads=nthreads())
   # assign indices, count observations and allocate a matrix
   Ψ = _alloc_lsq_matrix(configs, basis)
   # create the struct where everything is stored
   db = LsqDB(basis, configs, Ψ, dbpath)
   # parallel assembly of the LSQ matrix
   tfor_observations( configs,
      (n, okey, cfg, lck) -> safe_append!(db, lck, cfg, okey),
      msg = "Assemble LSQ blocks",
      verbose=verbose,
      maxnthreads=maxnthreads )
   # save to file
   if dbpath != ""
      verbose && @info("Writing db to disk...")
      try
         flush(db)
      catch
         @warn("""something went wrong trying to save the db to disk, but the data
               should be ok; if it is crucial to keep it, try to save manually.""")
      end
      verbose && @info("... done")
   else
      verbose && @info("db is not written to disk since `dbpath` is empty.")
   end
   return db
end

# function LsqDB_dist(basis::IPBasis,
#                      configs::AbstractVector{Dat})
#       # assign indices, count observations and allocate a matrix
#       Ψ = _alloc_lsq_matrix_dist(configs, basis, 4)
#       irow = 0
#       @sync @distributed for i = 1:nworkers()
#          # add energy to the lsq system
#          irow += 1
#          #y[irow] = wE * E / length(at)
#          #Ψ_part = localpart(Ψ) 
#          #@show size(Ψ_part)
#          #@show Ψ[irow, :] .= repeat([irow], 36)
#          Ψ_part = localpart(Ψ) 
#          vec(Ψ_part) .= (1:length(Ψ_part)) .+ 1000*myid()
#          #Ψ_part .= energy(basis, at.at) / length(at.at)
   
#          # add forces to the lsq system
#          # nf = 3*length(at.at)
#          # #y[(irow+1):(irow+nf)] = wF * mat(F)[:]
#          # Fb = forces(basis, at.at)
#          # for ib = 1:length(basis)
#          #    Ψ[(irow+1):(irow+nf), ib] = mat(Fb[ib])[:]
#          # end
#          # irow += nf
#       end
#       return Ψ
# end

function set_matrows!(d::Dat, okey::String, irows::Vector{Int})
   d.rows[okey] = irows
   return d
end

matrows(d::Dat, okey::String) = d.rows[okey]

function _alloc_lsq_matrix(configs, basis)
   # the indices associated with the basis are simply the indices within
   # the array - there is nothing else to do here
   #
   # loop through all observations and assign indices
   nrows = 0
   for (okey, d, _) in observations(configs)
      len = length(observation(d, okey))
      set_matrows!(d, okey, collect(nrows .+ (1:len)))
      nrows += len
   end
   # allocate and return the matrix
   return zeros(Float64, nrows, length(basis))
end

function _alloc_lsq_matrix_dist(configs, basis, nprocs)
   # the indices associated with the basis are simply the indices within
   # the array - there is nothing else to do here
   #
   # loop through all observations and assign indices
   nrows = 0
   for (okey, d, _) in observations(configs)
      len = length(observation(d, okey))
      set_matrows!(d, okey, collect(nrows .+ (1:len)))
      nrows += len
   end
   # allocate and return the matrix
   return dzeros((nrows,length(basis)), workers()[1:nprocs], [1,nprocs])
end


function safe_append!(db::LsqDB, db_lock, cfg, okey)
   # computing the lsq blocks ("rows") can be done in parallel,
   lsqrow = eval_obs(okey, basis(db), cfg)
   vec_lsqrow = vec_obs(okey, lsqrow)
   ### Something like this has to come here:
   if okey == "MU"
      tmp = zeros(3, length(vec_lsqrow))
      for i in 1:length(vec_lsqrow)
         tmp[:,i] = vec_lsqrow[i]
      end
      vec_lsqrow = tmp  
   end
   irows = matrows(cfg, okey)
   # but writing them into the DB must be done in a threadsafe way
   lock(db_lock)
   db.Ψ[irows, :] = vec_lsqrow
   unlock(db_lock)
   return nothing
end


# ------------------- Evaluating LSQ Blocks -----------------


# TODO: SOMETHING LIKE THIS IS STILL NEEDED...
# """
# evaluate one specific kind of datum, such as energy, forces, etc
# for one Dat across a large section of the basis; implicitly this will
# normally be done via `evaluate_many` or similar. The line
# ```
# vals = eval_obs(vDT, B, at)
# ```
# should likely be the bottleneck of `evallsq`
# """
# function _evallsq(vDT::Val,
#                   B::AbstractVector{<:AbstractCalculator},
#                   at::AbstractAtoms)
#    # vals will be a vector containing multiple evaluations
#    vals = eval_obs(vDT, B, at)
#    # vectorise the first so we know the length of the observation
#    vec1 = vec_obs(vDT, vals[1])
#    # create a multi-D array to reshape these into
#    A = Array{Float64}(undef, length(vec1), length(B))
#    A[:, 1] = vec1
#    for n = 2:length(B)
#       A[:, n] = vec_obs(vDT, vals[n])
#    end
#    return A
# end



# =========================== Convenience =================================


configtypes(db::LsqDB) = unique(configtype.(db.configs))

_nconfigs(db::LsqDB, ct::AbstractString) =
   length(find(configtype.(db.configs) .== ct))

# TODO: rewrite this!
function info(db::LsqDB)
   # config names, how many
   all_cts = configtype.(db.configs)
   unq_cts = unique(all_cts)
   configs_info = Dict{String, Int}()
   for ct in unq_cts
      configs_info[ct] = sum(all_cts .== ct)
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
   println(" Basis Group  | length | description ")
   for (i, B1) in enumerate(Bord)
      lenstr = replace(string(length(B1), pad=5), '0' => ' ')
      idxstr = replace(string(i, pad=2), '0' => ' ')
      desc = string(B1)
      desc = desc[1:min(50, length(desc))]
      println("          $idxstr  | $lenstr  | $desc" )
      #
   end
   println("======================================================")
end

end
