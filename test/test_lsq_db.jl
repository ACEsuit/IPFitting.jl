
using Test
using IPFitting, ProgressMeter, JuLIP, SHIPs
using JuLIP: decode_dict
using JuLIP.MLIPs: IPSuperBasis
Fit = IPFitting
DB = IPFitting.DB
Data = Fit.Data

function rand_data(sym, N, configtype="rand")
   cubic = (N == 1)
   at = bulk(sym, cubic=cubic) * N
   rattle!(at, 0.1)
   F = (N == 1) ? nothing : rand(3, length(at))
   return Dat(at, configtype; E = rand(), F = F) # , V = rand(3,3))
end

##
println("Double-Check (de-)dictionisation of basis: ")
basis1 = SHIPBasis(2, 15, 2.0, PolyTransform(2, 1.3), 2, 0.5, 3.0)
basis2 = SHIPBasis(3, 13, 2.0, PolyTransform(3, 1.0), 2, 0.5, 3.0)
B = IPSuperBasis(basis1, basis2)
println(@test decode_dict( Dict( basis1 ) ) == basis1)
println(@test decode_dict( Dict( basis2 ) ) == basis2)
println(@test decode_dict( Dict( B ) ) == B)

println("Double-Check (de-)dictionisation of Dat: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
println(@test Dat.(Dict.(data1)) == data1)
println(@test Dat.(Dict.(data2)) == data2)

##
println("Create a temporary database.")
tmpdir = mktempdir()
dbpath = joinpath(tmpdir, "temp")
data = [data1; data2]
db = nothing

try
   global db = DB.LsqDB(dbpath, B, data)
   println(@test true)
catch
   @info("...something went wrong...")
   println(@test false)
end

println("checking consistency of db")
println(@test DB.dbpath(db) == dbpath)
println(@test DB.kronfile(dbpath) == dbpath * "_kron.h5")
println(@test DB.infofile(dbpath) == dbpath * "_info.json")
println(@test db.basis == B)
println(@test db.configs == data)

println("re-load the database")
db1 = DB.LsqDB(dbpath)
println(@test DB.dbpath(db1) == dbpath)
println(@test db1.basis == db.basis)
println(@test db1.configs == db.configs)
println(@test db1.Ψ == db.Ψ)

##
println("re-load the database with mmap=true")
db2 = DB.LsqDB(dbpath, mmap=true)
println(@test DB.dbpath(db2) == dbpath)
println(@test db2.basis == db.basis)
println(@test db2.configs == db.configs)
println(@test db2.Ψ == db.Ψ)

##
println("Delete the temporary database")
rm(tmpdir; force=true, recursive=true)
