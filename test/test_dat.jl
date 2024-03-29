
using Test, JuLIP
using IPFitting
using IPFitting: Dat
using IPFitting.Data: configtype
using JuLIP: write_dict, read_dict
using JuLIP.Testing: print_tf 
##

@info("Testing Dat:")
println("  * Dat constructor")
println("  * vec, devec")
println("  * energy(dat), forces(dat), virial(dat)")
println("  * (de-)dictionisation, ==")

at = rattle!(bulk(:Fe) * 3, 0.1)
lj = LennardJones(r0 = rnn(:Fe)) * C2Shift(2.5*rnn(:Fe))
dat = Dat(at, "test"; E = energy(lj, at),
                      F = forces(lj, at),
                      V = virial(lj, at) )

print_tf(@test(configtype(dat) == "test"))
print_tf(@test(energy(dat) == energy(lj, at)))
print_tf(@test(forces(dat) == forces(lj, at)))
print_tf(@test(virial(dat) ≈ virial(lj, at)))
println() 

##

@info("Test (de-)serialisation of det")
print_tf(@test all(JuLIP.Testing.test_fio(dat)))
println() 