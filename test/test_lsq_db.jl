##
using Base.Test

using NBodyIPFitting
using NBodyIPs
using NBodyIPs.Polys: poly_basis
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
display(@test(DB._vec2arr(A) == A))
display(@test(DB._arr2vec(A) == A))
A = rand(3, 10)
B = [rand(3) for n = 1:10]
display(@test DB._vec2arr(DB._arr2vec(A)) == A)
display(@test DB._arr2vec(DB._vec2arr(B)) == B)
A = rand(3,4, 10)
B = [rand(3,4) for n = 1:10]
display(@test DB._vec2arr(DB._arr2vec(A)) == A)
display(@test DB._arr2vec(DB._vec2arr(B)) == B)

##
println("Check serialisation of basis: ")
basis1 = poly_basis(2, "exp( - 2 * (r/3-1))", "(:cos, 5.0, 7.0)", 10)
basis2 = poly_basis(3, "exp( - 2.5 * (r/3-1))", "(:cos2s, 2.0, 2.5, 4.0, 5.5)", 6)
display(@test Fit.Tools.decode.( Dict.( basis1 ) ) == basis1)
display(@test Fit.Tools.decode.( Dict.( basis2 ) ) == basis2)

println("Check serialisation of data: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
display(@test Dat.(Dict.(data1)) == data1)
display(@test Dat.(Dict.(data2)) == data2)


##
println("Create a temporary database.")
tmpdir = mktempdir()
db = DB.initdb(tmpdir, "db")
@test DB.dbdir(db) == tmpdir*"/db"
print("Check that initdb fails on existing database: ")
try
   initdb(tmpdir, "db")
   display(@test false)
catch
   display(@test true)
end

print("Add some basis functions: ")
try
   append!(db, basis1)
   append!(db, basis2)
   basis = DB.basis(db)
   display(@test true)
catch
   display(@test false)
end

print("Add some data to db: ")
try
   append!(db, data1)
   append!(db, data2)
   display(@test true)
catch
   display(@test false)
end

println("Reload the db")
db_dir = DB.dbdir(db)
db1 = DB.LsqDB(db_dir)
print("Test that the basis matches: ")
display(@test all(db.basis .== db1.basis))
print("Test that the data matches: ")
display(@test all(db.data .== db1.data))

println("Visually check that the files exist: " )
run(`ls $db_dir`)

println("Delete the temporary database")
rm(tmpdir; force=true, recursive=true)
