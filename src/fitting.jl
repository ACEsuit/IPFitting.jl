using JuLIP, ProgressMeter

export get_basis, regression, naive_sparsify,
       normalize_basis!, fiterrors

Base.norm(F::JVecsF) = norm(norm.(F))


"""
split off the inner assembly loop to
prepare for parallelising

* `d` : one data point
* `basis` : all basis functions
* `nforces` : randomly choose nforces
"""
function assemble_lsq_block(d, basis, nforces)
   at = Atoms(d)
   len = length(d)
   nforces = Int(min(nforces, len))
   # allocate observations (sub-) vector
   Y = zeros(1 + 3 * nforces)
   # allocate (sub-) matrix of basis functions
   Ψ = zeros(1 + 3 * nforces, length(basis))
   # ------- fill the data/observations vector -------------------
   Y[1] = energy(d)/len     # put in energy data
   if nforces > 0
      # extract the forces from the data:
      f = forces(d)
      If = rand(1:length(f), nforces)   # random subset of forces
      f_vec = mat(f[If])[:]             # convert it into a single long vector
      Y[2:end] = f_vec                  # put force data into rhs
   end
   # ------- fill the LSQ system, i.e. evaluate basis at data points -------
   for (ib, b) in enumerate(basis)
      Ψ[1, ib] = energy(b, at)/len
      # compute the forces
      if nforces > 0
         fb = forces(b, at)
         fb_vec = mat(fb[If])[:]
         Ψ[2:end, ib] = fb_vec
      end
   end
   # -------- what about the weight vector ------------
   return Ψ, Y
end

# TODO: parallelise!
function assemble_lsq(basis, data; verbose=true, nforces=0,
                      dt = verbose ? 0.5 : Inf)
   # generate many matrix blocks, one for each piece of data
   #  ==> this should be switched to pmap, or @parallel
   LSQ = @showprogress(dt, "assemble LSQ",
                  [assemble_lsq_block(d, basis, nforces) for d in data])
   # combine the local matrices into a big global matrix
   nY = sum(length(block[2]) for block in LSQ)
   Ψ = zeros(nY, length(basis))
   Y = zeros(nY)
   i0 = 0
   for id = 1:length(data)
      Ψi::Matrix{Float64}, Yi::Vector{Float64} = LSQ[id]
      rows = (i0+1):(i0+length(Yi))
      Ψ[rows, :] = Ψi
      Y[rows] = Yi
      i0 += length(Yi)
   end
   W = speye(length(Y))
   return Ψ, Y, I
end


function regression(basis, data;
                    verbose = true,
                    nforces=0, usestress=false,
                    stabstyle=:basis, stab=1e-3,
                    weights=:I)

   Ψ, Y, W = assemble_lsq(basis, data; verbose = verbose, nforces = nforces)
   if any(isnan, Ψ) || any(isnan, Y)
      error("discovered NaNs - something went wrong in the assembly")
   end

   @assert stabstyle == :basis

   # compute coefficients
   verbose && println("solve $(size(Ψ)) LSQ system using QR factorisation")
   Q, R = qr(Ψ)
   if W == I
      c = (R \ (Q' * Y)) ./ (1+stab)
   else
      A = Q' * (W * Q) + stab * eye(size(R, 1))
      b = Q' * (W * y)
      c = R \ (A \ b)
   end
   # check error on training set
   z = Ψ * c - Y
   rms = sqrt(dot(W * z, z) / length(Y))
   verbose && println("naive rms error on training set: ", rms)
   return c
end

# TODO:
#  - parallelise!
#  - combine rms and mae into one function
function fiterrors(V, data; verbose=true,
                   dt = verbose ? 0.5 : Inf)
   NE = 0
   NF = 0
   rmsE = 0.0
   rmsF = 0.0
   maeE = 0.0
   maeF = 0.0
   @showprogress dt "test" for d in data
      at, E, F = Atoms(d), energy(d), forces(d)
      # energy error
      Ex = energy(V, at)
      rmsE += (Ex - E)^2/length(at)^2
      maeE += abs(Ex-E) / length(at)
      NE += 1  # number of energies
      # force error
      Fx = forces(V, at)
      rmsF += sum( norm.(Fx - F).^2 )
      maeF += sum( norm.(Fx-F) )
      NF += length(Fx)   # number of forces
   end
   return sqrt(rmsE/NE), sqrt(rmsF/NF), maeE/NE, maeF / NF
end





"""
computes the maximum force over all configurations
in `data`
"""
function max_force(b, data)
   out = 0.0
   for d in data
      f = forces(b, Atoms(d))
      out = max(out, maximum(norm.(f)))
   end
   return out
end

function normalize_basis!(B, data)
   for b in B
      @assert (length(b.c) == 1)
      maxfrc = max_force(b, data)
      if maxfrc > 0
         b.c[1] /= maxfrc
      end
      if 0 < maxfrc < 1e-8
         warn("encountered a very small maxfrc = $maxfrc")
      end
   end
   return B
end

"""
remove a fraction `p` of normalised basis functions with the smallest
normalised coefficients.

NB: this returns complete *crap*.
"""
function naive_sparsify(B, c, data, p::AbstractFloat)
   # get normalisation constants for the basis functions
   nrmB = @showprogress "sparsify" [ force_norm(b, data) for b in B ]
   # normalised contributions
   cnrm = c .* nrmB
   @show nrmB
   # get the dominant contributions
   I = sortperm(abs.(cnrm))
   @show cnrm[I]
   # get the subset of indices to keep
   deleteat!(I, 1:floor(Int,length(B)*p))
   # return the sparse basis and corresponding coefficients
   return B[I], c[I]
end
