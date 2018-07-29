
##
using Base.Test
import NBodyIPFitting
import NBodyIPFitting.FIO

println("----------------------------")
println("Testing NBodyIPFitting.FIO")

fname1 = tempname() * ".h5"
fname2 = tempname() * ".h5"

D = Dict( "a" => rand(),
          "b" => rand(3,4,10),
          "c" => Dict( "d" => rand(10), "e" => "mystring"),
          "f" => nothing )

println("Save and Load Entire Dict")
FIO.save(fname1, D)
D1 = FIO.load(fname1)
println(@test D1 == D)

println("Load an element of D")
a = FIO.load(fname1, "a")
println(@test a == D["a"])

println("Save elements by name")
FIO.save(fname2, collect(D)...)
D2 = FIO.load(fname2)
println(@test D2 == D)

println("Read from multi-dim array")
bsub = FIO.load(fname1, ("b", (:,:,2:5)) )
println(@test bsub == D["b"][:,:,2:5])

println("Delete temporary files")
run(`rm $fname1 $fname2`)
println("----------------------------")
