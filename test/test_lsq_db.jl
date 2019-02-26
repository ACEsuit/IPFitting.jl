
using Test
using NBodyIPFitting, NBodyIPs, ProgressMeter, JuLIP
using NBodyIPs: blpolys
Fit = NBodyIPFitting
DB = NBodyIPFitting.DB
Data = Fit.Data

function rand_data(sym, N, configtype="rand")
   cubic = (N == 1)
   at = bulk(sym, cubic=cubic) * N
   rattle!(at, 0.1)
   F = (N == 1) ? nothing : rand(3, length(at))
   return Dat(at, configtype; E = rand(), F = F, V = rand(3,3))
end

##
println("Double-Check (de-)dictionisation of basis: ")
basis1 = blpolys(2, "exp( - 2 * (r/3-1))", "(:cos, 5.0, 7.0)", 10)
basis2 = blpolys(3, "exp( - 2.5 * (r/3-1))", "(:cos2s, 2.0, 2.5, 4.0, 5.5)", 6)
println(@test Fit.Tools.decode.( Dict.( basis1 ) ) == basis1)
println(@test Fit.Tools.decode.( Dict.( basis2 ) ) == basis2)

println("Double-Check (de-)dictionisation of Dat: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
println(@test Dat.(Dict.(data1)) == data1)
println(@test Dat.(Dict.(data2)) == data2)

##
println("Create a temporary database.")
tmpdir = mktempdir()
dbpath = joinpath(tmpdir, "temp")
basis = [basis1; basis2]
data = [data1; data2]
db = nothing
try
   global db = DB.LsqDB(dbpath, basis, data)
   println(@test true)
catch
   println("...something went wrong...")
   println(@test false)
end

println("checking consistency of db")
println(@test DB.dbpath(db) == dbpath)
println(@test DB.kronfile(dbpath) == dbpath * "_kron.h5")
println(@test DB.infofile(dbpath) == dbpath * "_info.json")
println(@test db.basis == basis)
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



# TODO: Some test like this will be needed when we switch to
#       loading the kron_groups on demand
# ## Confirm that the LSQ entries load correctly and
# db2 = DB.LsqDB(dbpath)
# println("Reload all dat files and check they are consistent with (data, basis)")
# @showprogress for n = 1:length(DB.data(db2))
#    lsq1 = DB.load_dat(db2, n)
#    lsq2 = Fit.Lsq.evallsq(DB.data(db2, 1), DB.basis(db2))
#    @test lsq1 == lsq2
# end
