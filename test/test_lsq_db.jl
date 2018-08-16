##
using Base.Test
using NBodyIPFitting, NBodyIPs, ProgressMeter
using NBodyIPs: blpolys
Fit = NBodyIPFitting
DB = NBodyIPFitting.DB
Data = Fit.Data

function rand_data(sym, N, config_type="rand")
   cubic = (N == 1)
   at = bulk(sym, cubic=cubic) * N
   rattle!(at, 0.1)
   F = (N == 1) ? nothing : rand(JVecF, length(at))
   return Dat(at, rand(), F, rand(JMatF), 1.0, config_type)
end

##
println("Testing _vec2arr and _arr2vec")
A = rand(10)
println(@test(DB._vec2arr(A) == A))
println(@test(DB._arr2vec(A) == A))
A = rand(3, 10)
B = [rand(3) for n = 1:10]
println(@test DB._vec2arr(DB._arr2vec(A)) == A)
println(@test DB._arr2vec(DB._vec2arr(B)) == B)
A = rand(3,4, 10)
B = [rand(3,4) for n = 1:10]
println(@test DB._vec2arr(DB._arr2vec(A)) == A)
println(@test DB._arr2vec(DB._vec2arr(B)) == B)

##
println("Check serialisation of basis: ")
basis1 = blpolys(2, "exp( - 2 * (r/3-1))", "(:cos, 5.0, 7.0)", 10)
basis2 = blpolys(3, "exp( - 2.5 * (r/3-1))", "(:cos2s, 2.0, 2.5, 4.0, 5.5)", 6)
println(@test Fit.Tools.decode.( Dict.( basis1 ) ) == basis1)
println(@test Fit.Tools.decode.( Dict.( basis2 ) ) == basis2)

println("Check serialisation of data: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
println(@test Dat.(Dict.(data1)) == data1)
println(@test Dat.(Dict.(data2)) == data2)


##
println("Create a temporary database.")
tmpdir = mktempdir()
db = DB.initdb(tmpdir, "db")
@test DB.dbdir(db) == tmpdir*"/db"
print("Check that initdb fails on existing database: ")
try
   initdb(tmpdir, "db")
   println(@test false)
catch
   println(@test true)
end

##
print("Add some basis functions: ")
# try
   append!(db, basis1)
   append!(db, basis2)
   basis = DB.basis(db)
   println(@test true)
# catch
#    println(@test false)
# end

print("Add some data to db: ")
# try
append!(db, data1)
append!(db, data2)
println(@test true)
# catch
#    println(@test false)
# end

println("Reload the db")
db_dir = DB.dbdir(db)
db1 = DB.LsqDB(db_dir)
print("Test that the basis matches: ")
println(@test all(db.basis .== db1.basis))
print("Test that the data matches: ")
println(@test all(db.data .== db1.data))

# ##
# println("Add additional basis functions")
# basis3 = blpolys(3, "(3/r)^5", "(:cos2s, 2.0, 2.5, 3.5, 5.5)", 5)
# append!(db, basis3)
#
# ## Confirm that the LSQ entries load correctly and
# db2 = DB.LsqDB(db_dir)
# println("Reload all dat files and check they are consistent with (data, basis)")
# @showprogress for n = 1:length(DB.data(db2))
#    lsq1 = DB.load_dat(db2, n)
#    lsq2 = Fit.Lsq.evallsq(DB.data(db2, 1), DB.basis(db2))
#    @test lsq1 == lsq2
# end

# lsq2 = Fit.Lsq.evallsq(DB.data(db2, 1), DB.basis(db2))
# db2.data[1]
# lsq1 = DB.load_dat(db1, 1)

##
println("Delete the temporary database")
rm(tmpdir; force=true, recursive=true)
