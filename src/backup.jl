


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


      # TODO: this should be generalised to be able to hook into
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




# include("poly_regularise.jl")


# # ----------------- some simplified access functions ------------------
#
# evaluate(V::NBodyFunction{2}, r::Number) = evaluate(V, SVector(r))
#
# evaluate_d(V::NBodyFunction{2}, r::Number) = evaluate_d(V, SVector(r))[1]
#
# evaluate_dd(V::NBodyFunction{2}, r::Number) =
#       ((@D V(r+1e-5)) - (@D V(r-1e-5))) / 1e-5
#
# evaluate(V::NBodyFunction{3}, r1::Number, r2::Number, r3::Number) =
#       evaluate(V, SVector(r1, r2, r3))



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
