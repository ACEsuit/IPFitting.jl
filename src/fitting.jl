
using JuLIP

export get_basis, regression

function get_basis(ord, dict, sym, rcut;
                   degree=:default, kwargs...)
   dim = (ord * (ord-1)) ÷ 2
   if degree == :default
      _, fs, dfs = psym_polys(dim, dict, sym; kwargs...)
   elseif degree == :total
      _, fs, dfs = psym_polys_tot(dim, dict, sym; kwargs...)
   else
      error("unkown degree type")
   end
   return NBody.(ord, fs, dfs, rcut)
end


function regression(basis, data; verbose = true)
   A = zeros(length(data), length(basis))
   F = zeros(length(data))
   println("assemble system")
   lenat = 0
   for (id, d) in enumerate(data)
      at = d[1]::Atoms{Float64, Int}
      lenat = max(lenat, length(at))
      F[id] = d[2]::Float64
      for (ib, b) in enumerate(basis)
         print(".")
         A[id, ib] = b(at)
      end
   end
   # compute coefficients
   println("solve lsq")
   Q, R = qr(A)
   c = R \ (Q' * F)
   # check error on training set
   println("rms error on training set: ",
           norm(A * c - F) / sqrt(length(data)) / sqrt(lenat) )
   return c
end
