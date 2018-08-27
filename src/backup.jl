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


function regularise(Ψ, Y, P::Matrix)
   @assert size(Ψ,2) == size(P,2)
   return vcat(Ψ, P), vcat(Y, zeros(size(P,1)))
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


"""

## Keyword Arguments:

* weights: either `nothing` or a tuple of `Pair`s, i.e.,
```
weights = (:E => 100.0, :F => 1.0, :V => 0.01)
```
Here `:E` stand for energy, `:F` for forces and `:V` for virial .

* config_weights: a tuple of string, value pairs, e.g.,
```
config_weights = ("solid" => 10.0, "liquid" => 0.1)
```
this adjusts the weights on individual configurations from these categories
if no weight is provided then the weight provided with the is used.
Note in particular that `config_weights` takes precedence of Dat.w!
If a weight 0.0 is used, then those configurations are removed from the LSQ
system.
"""


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
