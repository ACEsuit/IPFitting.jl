
using JuLIP, ProgressMeter

export get_basis, regression, rms

Base.norm(F::JVecsF) = norm(norm.(F))

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

function assemble_lsq(basis, data; verbose=true, nforces=0)
   A = zeros(length(data) * (1+3*nforces), length(basis))
   F = zeros(length(data) * (1+3*nforces))
   lenat = 0
   if verbose
      pm = Progress(length(data), desc="assemble LSQ system")
   end
   for (id, d) in enumerate(data)
      i0 = (id-1) * (1+3*nforces) + 1
      at = d[1]
      lenat = max(lenat, length(at))
      F[i0] = d[2]::Float64     # put in energy data
      if nforces > 0
         # extract the force from the data:
         f = d[3]::JVecsF    # a vector of short vectors
         If = rand(1:length(f), nforces)   # random subset of forces
         f_vec = mat(f[If])[:]   # convert it into a single long vector
         @assert !any(isnan.(f_vec))
         F[(i0+1):(i0+3*nforces)] = f_vec    # put force data into rhs
      end
      for (ib, b) in enumerate(basis)
         A[i0, ib] = b(at)
         # compute the forces
         if nforces > 0
            fb = - @D b(at)
            fb_vec = mat(fb[If])[:]
            @assert !any(isnan.(fb_vec))
            A[(i0+1):(i0+3*nforces), ib] = fb_vec
         end
      end
      if verbose
         next!(pm)
      end
   end
   return A, F, lenat
end

function regression(basis, data; verbose = true, nforces=0)
   A, F, lenat = assemble_lsq(basis, data;
                     verbose = verbose, nforces = nforces)
   @assert !any(isnan.(A))
   # compute coefficients
   verbose && println("solve $(size(A)) LSQ system using QR factorisation")
   Q, R = qr(A)
   c = R \ (Q' * F)
   # check error on training set
   verbose && println("rms error on training set: ",
                       norm(A * c - F) / sqrt(length(data)) / sqrt(lenat) )
   return c
end


function rms(V, data)
   NE = 0
   NF = 0
   errE = 0.0
   errF = 0.0
   for (at, E, F) in data
      # energy error
      Ex = energy(V, at)
      errE += (Ex - E)^2
      NE += 1
      # force error
      Fx = forces(V, at)
      errF += sum( norm.(Fx - F)^2 )
      NF += length(Fx)
   end
   return sqrt(errE/NE), sqrt(errF/NF)
end


function mae(V, data)
   NE = 0
   NF = 0
   errE = 0.0
   errF = 0.0
   for (at, E, F) in data
      # energy error
      Ex = energy(V, at)
      errE += abs(Ex - E)
      NE += 1
      # force error
      Fx = forces(V, at)
      errF += sum( abs.(Fx - F) )
      NF += length(Fx)
   end
   return errE / NE, errF / NF
end
