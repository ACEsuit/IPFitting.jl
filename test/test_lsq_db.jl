
using Test
using IPFitting, ProgressMeter, JuLIP, ACE1
using JuLIP: read_dict
using JuLIP.Testing: test_fio, print_tf
using JuLIP.MLIPs: IPSuperBasis
Fit = IPFitting
DB = IPFitting.DB
Data = Fit.Data

##
   
function rand_data(sym, N, configtype="rand")
   cubic = (N == 1)
   at = bulk(sym, cubic=cubic) * N
   rattle!(at, 0.1)
   F = (N == 1) ? nothing : rand(3, length(at))
   return Dat(at, configtype; E = rand(), F = F) # , V = rand(3,3))
end


##
@info("Double-Check (de-)dictionisation of basis: ")
basis1 = rpi_basis(species = :Ti, N = 2, maxdeg = 8)
basis2 = rpi_basis(species = :Ti, N = 3, maxdeg = 6)
B = IPSuperBasis(basis1, basis2)
print_tf(@test all(test_fio(basis1)))
print_tf(@test all(test_fio(basis2)))
print_tf(@test all(test_fio(B)))
println() 

@info("Double-Check (de-)dictionisation of Dat: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
print_tf(@test read_dict.(write_dict.(data1)) == data1)
print_tf(@test read_dict.(write_dict.(data2)) == data2)
println() 

##
@info("Create a temporary database.")
dbpath = tempname()
data = [data1; data2]
db = nothing

try
   global db = DB.LsqDB(dbpath, B, data)
   print_tf(@test true)
catch
   @info("...something went wrong...")
   print_tf(@test false)
end
println() 
@info("checking consistency of db")
print_tf(@test DB.dbpath(db) == dbpath)
print_tf(@test DB.kronfile(dbpath) == dbpath * "_kron.h5")
print_tf(@test DB.infofile(dbpath) == dbpath * "_info.json")
print_tf(@test db.basis == B)
print_tf(@test db.configs == data)
println() 

@info("re-load the database")
db1 = DB.LsqDB(dbpath)
print_tf(@test DB.dbpath(db1) == dbpath)
print_tf(@test db1.basis == db.basis)
print_tf(@test db1.configs == db.configs)
print_tf(@test db1.Ψ == db.Ψ)
println() 

##
@info("re-load the database with mmap=true")
db2 = DB.LsqDB(dbpath, mmap=true)
print_tf(@test DB.dbpath(db2) == dbpath)
print_tf(@test db2.basis == db.basis)
print_tf(@test db2.configs == db.configs)
print_tf(@test db2.Ψ == db.Ψ)
println() 