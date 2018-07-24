
module DB

using JuLIP, ProgressMeter, Base.Threads

using NBodyIPFitting: Dat, LsqDB

import Base: append!, push!

export get_basis,
       observations, get_lsq_system,
       regularise, table, load_lsq, load_ip,
       config_types,
       del_data!

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
   dbdir = jointpath(basedir, dbname)
   @assert !isdir(dbdir)
   mkdir(dbdir)
   save(joinpath(dbdir, "data.jld2"), "data", Dat{Float64}[])
   save(joinpath(dbdir, "basis.jld2"), "basis", AbstractCalculator[])
   return LsqDB(dbdir)
end

function LsqDB(dirname::AbstractString)
   if !isdir(dirname)
      error("""`dirname` is not a directory. If you want to create a new LsqDB
               please use `initdb`.""")
   end

end

push!(db::LsqDB, d::Dat) = append!(db::LsqDB, [d])

function append!(db::LsqDB, d::AbstractVector{TD}) where {TD <: Dat}
end

push!(db::LsqDB, b::AbstractCalculator) = append!(db, [b])

function append!(db::LsqDB, b::AbstractVector{TB}) where {TB <: AbstractCalculator}
end



end
