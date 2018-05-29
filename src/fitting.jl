using JuLIP, ProgressMeter

export get_basis, regression, naive_sparsify,
       normalize_basis!, fiterrors, scatter_data

Base.norm(F::JVecsF) = norm(norm.(F))

# components of the stress (up to symmetry)
const _IS = SVector(1,2,3,5,6,9)




function assemble_lsq_block(d, Bord, Iord, nforces,
                            w_E, w_F, w_S)
   len = length(d)
   nforces = Int(min(nforces, len))
   # ------- fill the data/observations vector -------------------
   Y = Float64[]
   # energy
   push!(Y, w_E * energy(d))
   # forces
   if forces(d) != nothing
      f = forces(d)
      # If = rand(1:length(f), nforces)   # random subset of forces
      # f_vec = mat(f[If])[:]             # convert it into a single long vector
      f_vec = mat(f)[:]
      append!(Y, f_vec)                 # put force data into rhs
   end
   # stress / virial
   if virial(d) != nothing
      S = virial(d)
      append!(Y, w_S * S[_IS])
   end

   # ------- fill the LSQ system, i.e. evaluate basis at data points -------
   at = Atoms(d)
   # allocate (sub-) matrix of basis functions
   Ψ = zeros(length(Y), sum(length.(Bord)))

   # energies
   i0 = 0
   for n = 1:length(Bord)
      Es = energy(Bord[n], at)
      Ψ[i0+1, Iord[n]] = w_E * Es
   end
   i0 += 1

   # forces
   if forces(d) != nothing
      for n = 1:length(Bord)
         Fs = forces(Bord[n], at)
         for j = 1:length(Fs)
            # fb_vec = mat(Fs[j][If])[:]
            fb_vec = mat(Fs[j])[:]
            Ψ[(i0+1):(i0+length(fb_vec)), Iord[n][j]] = fb_vec
         end
      end
   end
   i0 += 3 * nforces

   # stresses
   if virial(d) != nothing
      for n = 1:length(Bord)
         Ss = virial(Bord[n], at)
         for j = 1:length(Ss)
            Svec = w_S * Ss[j][_IS]
            Ψ[(i0+1):(i0+length(_IS)), Iord[n][j]] = Svec
         end
      end
   end

   # -------- what about the weight vector ------------
   return d.w * Ψ, d.w * Y
end


function split_basis(basis)
   # get the types of the individual basis elements
   tps = typeof.(basis)
   Iord = Vector{Int}[]
   Bord = Any[]
   for tp in unique(tps)
      # find which elements of basis have type `tp`
      I = find( tp .== tps )
      push!(Iord, I)
      push!(Bord, [ b for b in basis[I]] )
   end
   return Bord, Iord
end



# TODO: parallelise!
function assemble_lsq(basis, data; verbose=true, nforces=Inf,
                      dt = verbose ? 0.5 : Inf,
                      w_E = 200.0, w_F = 1.0,
                      w_S = 0.3 )
   # sort basis set into body-orders, and possibly different
   # types within the same body-order (e.g. for matching!)
   Bord, Iord = split_basis(basis)

   # generate many matrix blocks, one for each piece of data
   #  ==> this should be switched to pmap, or @parallel
   LSQ = @showprogress(dt, "assemble LSQ",
               [assemble_lsq_block(d, Bord, Iord, nforces, w_E, w_F, w_S)
                for d in data])
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
                    stabstyle=:basis, cstab=1e-3,
                    weights=:I,
                    regulariser = nothing,
                    w_E = 30.0, w_F = 1.0, w_S = 0.3)

   Ψ, Y, W = assemble_lsq(basis, data;
               verbose = verbose, nforces = nforces,
               w_E = w_E, w_F = w_F, w_S = w_S)

   if any(isnan, Ψ) || any(isnan, Y)
      error("discovered NaNs - something went wrong in the assembly")
   end

   @assert stabstyle == :basis

   # compute coefficients
   verbose && println("solve $(size(Ψ)) LSQ system using QR factorisation")
   Q, R = qr(Ψ)
   @show cond(R)

   if weights == :E
      ndat = length(Y) ÷ length(data)
      w = zeros(ndat, length(data))
      E0 = minimum(energy(d) for d in data)
      for (n, d) in enumerate(data)
         w[:, n] = 1.0 / ( energy(d) / length(d) - E0 + 1.0 )^2
         W = Diagonal(w[:])
      end
   end

   if W == I && regulariser == nothing
      c = (R \ (Q' * Y)) ./ (1+cstab)
   elseif regulariser == nothing
      A = Q' * (W * Q) + cstab * eye(size(R, 1))
      b = Q' * (W * Y)
      c = R \ (A \ b)
   else
      @assert W == I
      A = (1 + cstab) * R' * R + regulariser
      b = R' * Q' * Y
      c = A \ b
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
   @showprogress dt "fiterrors" for d in data
      at, E, F = Atoms(d), energy(d), forces(d)
      # energy error
      Ex = energy(V, at)
      rmsE += (Ex - E)^2 / length(at)^2
      maeE += abs(Ex-E) / length(at)
      NE += 1  # number of energies
      # force error
      Fx = forces(V, at)
      rmsF += sum( norm.(Fx - F).^2 )
      maeF += sum( norm.(Fx-F) )
      NF += 3*length(Fx)   # number of forces
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



function scatter_data(IP, data)
   E_data = Float64[]
   E_fit = Float64[]
   F_data = Float64[]
   F_fit = Float64[]
   for d in data
      at = Atoms(d)
      len = length(at)
      push!(E_data, energy(d) / len)
      append!(F_data, mat(forces(d))[:])
      push!(E_fit, energy(IP, at) / len)
      append!(F_fit, mat(forces(IP, at))[:])
   end
   return E_data, E_fit, F_data, F_fit
end
