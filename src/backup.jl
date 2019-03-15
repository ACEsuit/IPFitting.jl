function hess_weights_hook!(w, d::Dat)
   at = Atoms(d)
   if energy(d) == nothing || forces(d) == nothing
      warn("""hess_weights_hook!: a training configuration does not contain
              energy and forces => ignore""")
     return w
  end
   # compute R vectors
   X = positions(at)
   h = 0.01 # norm(X[1])
   # if h < 1e-5 || h > 0.02
   #    warn("unexpected location of X[1] in hess_weights_hook => ignore")
   #    @show X[1], h
   #    return w
   # end

   # give this energy a lot of weight to make sure we match the
   # ground state (which we assume this is)
   w[1] = 0.0

   # now fix the scaling for the force weights
   # X[1] *= 0
   R = [ JuLIP.project_min(at, x - X[1])  for x in X ]
   r = norm.(R)
   r3 = (ones(3) * r')[:]
   If = 2:(1+3*length(R))
   w[If] .= w[If] .* (r3.^7) / h
   w[2] = 0.0
   return w
end



using NBodyIPs.Data: config_type

function Base.info(lsq::LsqSys)
   println(repeat("=", 60))
   println(" LsqSys Summary")
   println(repeat("-", 60))
   println("      #configs: $(length(lsq.data))")
   println("    #basisfcns: $(length(lsq.basis))")
   println("  config_types: ",
         prod(s*", " for s in config_types(lsq)))

   Bord, _ = split_basis(lsq.basis)
   println(" #basis groups: $(length(Bord))")
   println(repeat("-", 60))

   for (n, B) in enumerate(Bord)
      println("   Group $n:")
      info(B; indent = 6)
   end
   println(repeat("=", 60))
end


function Base.append!(lsq::NBodyIPs.LsqSys, data::Vector{TD}) where TD <: NBodyIPs.Data.Dat
   return append!(lsq, kron(data, lsq.basis))
end

function Base.append!(lsq::NBodyIPs.LsqSys, lsq2::NBodyIPs.LsqSys)
   lsq.data = [lsq.data; lsq2.data]
   lsq.Ψ = [lsq.Ψ; lsq2.Ψ]
   return lsq
end




function analyse_subbasis(lsq, order, degrees, Ibasis)
   not_nothing = length(find( [order, degrees, Ibasis] .!= nothing ))
   if not_nothing > 1
      @show order, degrees, Ibasis
      error("at most one of the keyword arguments `order`, `degrees`, `Ibasis` maybe provided")
   end

   # if the user has given no information, then we just take the entire basis
   if not_nothing == 0
      order = Inf
   end

   # if basis is chosen by maximum body-order ...
   if order != nothing
      Ibasis = find(1 .< bodyorder.(lsq.basis) .<= order) |> sort

   # if basis is chosen by maximum degrees body-order ...
   # degrees = array of
   elseif degrees != nothing
      # check that the constant is excluded
      @assert degrees[1] == 0
      Ibasis = Int[]
      for (deg, I) in zip(degrees, lsq.Iord)
         if deg == 0
            continue
         end
         # loop through all basis functions in the current group
         for i in I
            # and add all those to Ibasis which have degree <= deg
            if degree(lsq.basis[i]) <= deg
               push!(Ibasis, i)
            end
         end
      end

   # of indeed if the basis is just constructed through individual indices
   else
      Ibasis = Int[i for i in Ibasis]
   end

   return Ibasis
end



* `order`: integer specifying up to which body-order to include the basis
in the fit. (default: all basis functions are included)

normalise_E = true, normalise_V = true,
hooks = Dict("hess" => hess_weights_hook!)) =


      #       this should be generalised to be able to hook into
      #       other kinds of weight modifications
      idx_end = idx
      if w != 0 && length(config_type(d)) >= 4 && config_type(d)[1:4] == "hess"
         if haskey(hooks, "hess")
            w = hooks["hess"](W[idx_init:idx_end], d)
            W[idx_init:idx_end] = w
         end
      end
   end




# ------- distribution functions on the invariants --------------

function rdf(at::Atoms, rcut, transform=identity)
   if transform != identity
      return idf(2, at, rcut, transform)[1][1]
   end
   nlist = neighbourlist(at, rcut)
   return nlist.r
end

"""
`idf(N::Integer, at, rcut, transform)`

invariants distribution function => accumulates all the invariants values
arising during an N-body assembly over at.
"""
idf(N::Integer, at, rcut, transform) =
      _idf(Val(N), Val(bo2edges(N)), at, rcut, transform)

function _idf(valN::Val{N}, valM::Val{M}, at::Atoms{T}, rcut::T,
              transform) where {N, M, T}
   # compute invariants vectors to learn how many there are
   x = rand(SVector{M,T})
   I1, I2 = invariants(x)

   I1acc = [ T[] for n = 1:length(I1) ]
   I2acc = [ T[] for n = 1:length(I2) ]
   Iacc = (I1acc, I2acc)

   for (i, j, r, R) in sites(at, rcut)
      eval_site_nbody!(valN, R, rcut,
                  (Iacc, s, _1, _2, _3) -> idf_accumulator(Iacc, transform.(s)),
                  Iacc, nothing)
   end

   return I1acc, I2acc
end

function idf_accumulator(Iacc, x)
   I1, I2 = invariants(x)
   for n = 1:length(I1)
      push!(Iacc[1][n], I1[n])
   end
   for n = 1:length(I2)
      push!(Iacc[2][n], I2[n])
   end
   return Iacc
end




# ---------------------- Regularisation of EnvIPs --------------------------

function regularise(N::Integer, P::Integer, B::Vector, r0, r1; kwargs...)
   I = find(isa.(B, EnvBL{N, P}))
   Br = [b.Vr for b in B[I]]
   Ψr = Polys.regularise(N, Br, r0, r1; kwargs...)
   Ψ = zeros(size(Ψr, 1), length(B))
   Ψ[:, I] = Ψr
   return Ψ
end

function regularise(N::Integer, B::Vector, r0, r1; kwargs...)
   # find all EnvBL basis functions with body-order N
   I = find(isa.(B, EnvBL{N}))
   max_P = maximum( b.t for b in B[I] )
   # each polynomial degree p => separate V_Np, hence
   # compute multiple regularisations
   Ψ = zeros(0, length(B))
   for p = 0:max_P
      Ψp = regularise(N, p, B, r0, r1; kwargs...)
      Ψ = [Ψ; Ψp]
   end
   return Ψ
end


# # =============== Experimental:
# #   evaluate NBodyIP
#
# (V::NBodyIP)(args...) = evaluate(V, args...)
#
# evaluate(V::NBodyIP, r::Number) = evaluate(V::NBodyIP, SVector(r))
#
# evaluate(V::NBodyIP, r1::T, r2::T, r3::T) where {T <: Number} =
#       evaluate(V::NBodyIP, SVector(r1, r2, r3))
#
# function evaluate(V::NBodyIP, r::SVector{N, T}) where {N, T}
#    v = zero(T)
#    for Vn in V.components
#       if bo2edges(bodyorder(Vn)) == N
#          v += Vn(r)
#       end
#    end
#    return v
# end

# dim(V::NBPoly{N,M}) where {N, M} = M-1
