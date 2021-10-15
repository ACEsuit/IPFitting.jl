
using Test
using IPFitting, ProgressMeter, JuLIP, ACE, ACEatoms
using ACE: SymmetricBasis, SimpleSparseBasis
using JuLIP: read_dict
using JuLIP.Testing: test_fio
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
println("Double-Check (de-)dictionisation of basis: ")
maxdeg1 = 8
ord1 = 2
Bsel1 = SimpleSparseBasis(ord1, maxdeg1)
B1p1 = ACEatoms.ZμRnYlm_1pbasis(; species = [:Ti], maxdeg = maxdeg1, Bsel = Bsel1, 
                                 rin = 1.2, rcut = 5.0)
ACE.init1pspec!(B1p1, Bsel1)
basis1 = SymmetricBasis(ACE.Invariant(), B1p1, Bsel1)

maxdeg2 = 6
ord2 = 3
Bsel2 = SimpleSparseBasis(ord2, maxdeg2)
B1p2 = ACEatoms.ZμRnYlm_1pbasis(; species = [:Ti], maxdeg = maxdeg2, Bsel = Bsel2, 
                                 rin = 1.2, rcut = 5.0)
ACE.init1pspec!(B1p2, Bsel2)
basis2 = SymmetricBasis(ACE.Invariant(), B1p2, Bsel2)
B = IPSuperBasis([basis1, basis2])

println(@test all(test_fio(basis1)))
println(@test all(test_fio(basis2)))
#println(@test all(test_fio(B)))  # This doesn't work !!

println("Double-Check (de-)dictionisation of Dat: ")
data1 = [ rand_data(:Ti, 3, "md") for n = 1:10 ]
data2 = [ rand_data(:Ti, 1, "cell") for n = 1:10 ]
println(@test read_dict.(write_dict.(data1)) == data1)
println(@test read_dict.(write_dict.(data2)) == data2)

##
println("Create a temporary database.")
dbpath = tempname()
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
