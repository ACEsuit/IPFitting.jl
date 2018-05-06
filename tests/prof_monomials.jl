using BenchmarkTools, StaticArrays

include("fastpolys.jl")
using FastPolys

const I1 = (1,1,1,2,1,1,1,2,2,3,4,3,4,2,3,4,5,5,6,7)
const I2 = (2,2,3,3,5,5,6,5,5,6,7,6,7,8,8,9,6,8,8,9)
const I3 = (4,3,4,4,7,6,7,9,8,10,10,8,9,9,10,10,7,9,10,10)
const I4 = (4,4,4 ,4,4,4,7,7,7 ,8,8,8 ,9 ,8 ,9 ,9 ,9 ,9 ,6,6)
const I5 = (8,9,10,6,7,5,8,9,10,9,9,10,10,10,10,10,10,10,7,8)
const II = [I1,I2,I3,I4,I5]

# const I12 = Val((I1, I2))
# const I123 = Val((I1, I2, I3))
# const I1234 = Val((I1, I2, I3, I4))
# const I12345 = Val((I1, I2, I3, I4, I5))

x = [(@SVector rand(10)) for i = 1:5]
dx = [(@SVector rand(10)) for i = 1:5]

for n = 1:5
   println("Order = $n:")
   X = (x[1:n]...)
   dX = (dx[1:n]...)
   P = Val((II[1:n]...))
   print("    fpoly: "); @btime fpoly($X, $P)
   print("  fpoly_d: "); @btime fpoly_d($X, $dX, $P)
end
