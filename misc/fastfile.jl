
import JLD2, FileIO
using Serialization, Test, HDF5

for N in [1_000, 2_000, 3_000, 4_000]

A = rand(N,200_000)
sizeof(A) / 1e9

@info "N = $N"
@info "sizeof(A) = $(round(sizeof(A) / 1e9, 2))GB"
@info "Testing JLD2"
@time FileIO.save("temp.jld2", "A", A)
@time B = FileIO.load("temp.jld2", "A")
println(@test B == A)
run(`rm temp.jld2`)

@info "Testing serialize/deserialize"
function mysave(A, fname)
   f = open(fname, "w")
   serialize(f, A)
   close(f)
end
function myload(fname)
   f = open(fname, "r")
   A = deserialize(f)
   close(f)
   return A
end
@time mysave(A, "temp.dat")
@time B = myload("temp.dat")
run(`rm temp.dat`)
println(@test A == B)


@info "Testing HDF5"
function savehdf5(A, fname)
   h5open(fname, "w") do fid
      fid["A"] = A
   end
end
function loadhdf5(fname)
   fid = h5open(fname, "r")
   A = read(fid["A"])
   close(fid)
   return A
end
@time savehdf5(A, "temp.h5")
@time B = loadhdf5("temp.h5")
run(`rm temp.h5`)
println(@test A == B)

end
