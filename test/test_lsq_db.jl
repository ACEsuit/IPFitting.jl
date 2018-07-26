##
using Base.Test

using NBodyIPFitting
using NBodyIPs
using NBodyIPs.Polys: poly_basis
Fit = NBodyIPFitting
DB = NBodyIPFitting.DB

##
println("Testing _vec2arr and _arr2vec")
A = rand(10)
@test DB._vec2arr(A) == A
@test DB._arr2vec(A) == A
A = rand(3, 10)
B = [rand(3) for n = 1:10]
@test DB._vec2arr(DB._arr2vec(A)) == A
@test DB._arr2vec(DB._vec2arr(B)) == B
A = rand(3,4, 10)
B = [rand(3,4) for n = 1:10]
@test DB._vec2arr(DB._arr2vec(A)) == A
@test DB._arr2vec(DB._vec2arr(B)) == B

##
println("Create a temporary database")
tmpdir = mktempdir()
db = DB.initdb(tmpdir, "db")
@test DB.dbdir(db) == tmpdir*"/db"
println("Check that initdb fails on existing database.")
try
   initdb(tmpdir, "db")
   display(@test false)
catch
   display(@test true)
end

println("Add some basis functions")
try
   append!(db, poly_basis(2, "exp( - 2 * (r/3-1))", "(:cos, 5.0, 7.0)", 10))
   append!(db, poly_basis(3, "exp( - 2.5 * (r/3-1))", "(:cos2s, 2.0, 2.5, 4.0, 5.5)", 6))
   basis = DB.basis(db)
   display(@test true)
catch
   display(@test false)
end

println("Reload the db")
db_dir = DB.dbdir(db)
db1 = DB.LsqDB(db_dir)

# COMPARE OLD DB AND NEW DB => NEED COMPARISON OF THE basis functions

rm(tmpdir; force=true, recursive=true)
