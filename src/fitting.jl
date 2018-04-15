using JuLIP, ProgressMeter

export get_basis, regression, rms, mae

Base.norm(F::JVecsF) = norm(norm.(F))

# split off the inner assembly loop to
# prepare for parallelising
function assemble_lsq_block(d, basis, nforces)
   F = zeros(1 + 3 * nforces)
   A = zeros(1 + 3 * nforces, length(basis))
   at = d[1]::Atoms
   # ---- fill the data vector -------------------
   F[1] = d[2]::Float64     # put in energy data
   if nforces > 0
      # extract the forces from the data:
      f = d[3]::JVecsF                  # a vector of short vectors
      If = rand(1:length(f), nforces)   # random subset of forces
      f_vec = mat(f[If])[:]             # convert it into a single long vector
      F[2:end] = f_vec                  # put force data into rhs
   end
   # ------- fill the LSQ system, i.e. evaluate basis at data points -------
   for (ib, b) in enumerate(basis)
      A[1, ib] = energy(b, at)
      # compute the forces
      if nforces > 0
         fb = forces(b, at)
         fb_vec = mat(fb[If])[:]
         A[2:end, ib] = fb_vec
      end
   end
   return (A, F, length(at))
end

function assemble_lsq(basis, data; verbose=true, nforces=0)
   A = zeros(length(data) * (1+3*nforces), length(basis))
   F = zeros(length(data) * (1+3*nforces))
   lenat = 0
   # generate many matrix blocks, one for each piece of data
   #  ==> this should be switched to pmap!
   if verbose
      LSQ = @showprogress [assemble_lsq_block(d, basis, nforces) for d in data]
   else
      LSQ = [assemble_lsq_block(d, basis, nforces) for d in data]
   end
   # lsq_block = d -> assemble_lsq_block(d, basis, nforces)
   # LSQ = pmap(lsq_block, data, distributed=false)
   # combine the local matrices into a big global matrix
   for id = 1:length(data)
      i0 = (id-1) * (1+3*nforces) + 1
      rows = i0:(i0+3*nforces)
      A_::Matrix{Float64}, F_::Vector{Float64}, lenat_::Int = LSQ[id]
      lenat = max(lenat, lenat_)
      A[rows, :] = A_
      F[rows] = F_
   end
   return A, F, lenat
end

# TODO: parallelise!
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
                       norm(A * c - F) / sqrt(length(F)) )
   return c
end

# TODO:
#  - parallelise!
#  - combine rms and mae into one function 
function rms(V, data)
   NE = 0
   NF = 0
   errE = 0.0
   errF = 0.0
   @showprogress for n = 1:length(data)
      at, E, F = data[n]
      # energy error
      Ex = energy(V, at)
      errE += (Ex - E)^2
      NE += length(at)   # number of site energies
      # force error
      Fx = forces(V, at)
      errF += sum( norm.(Fx - F).^2 )
      NF += length(Fx)   # number of forces
   end
   return sqrt(errE/NE), sqrt(errF/NF)
end


# TODO: parallelise!
function mae(V, data)
   NE = 0
   NF = 0
   errE = 0.0
   errF = 0.0
   @showprogress for n = 1:length(data)
      at, E, F = data[n]
      # energy error
      Ex = energy(V, at)
      errE += abs(Ex - E)
      NE += length(at)   # number of site energies
      # force error
      Fx = forces(V, at)
      errF += sum( norm.(Fx - F) )
      NF += length(Fx)   # number of forces
   end
   return errE/NE, errF/NF
end
