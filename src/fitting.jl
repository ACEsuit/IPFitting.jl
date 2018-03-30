
using JuLIP

export get_basis, regression

function get_basis(ord, dict, sym, rcut;
                   degree=:default, kwargs...)
   if degree == :default
      _, fs, dfs = psym_polys(ord, dict, sym; kwargs...)
   elseif degree == :total
      _, fs, dfs = psym_polys_tot(ord, dict, sym; kwargs...)
   else
      error("unkown degree type")
   end
   return NBody.(ord, fs, dfs, rcut)
end


function regression(basis, data)
   A = zeros(length(data), length(basis))
   F = zeros(length(data))
   for (id, d) in enumerate(data)
      at = d[1]::Atoms{Float64, Int}
      F[id] = d[2]::Float64
      A[id, :] = basis.(at)
   end
   # compute coefficients
   Q, R = qr(A)
   return R \ (Q' * F)
end
