using JuLIP, ProgressMeter
using NBodyIPs.Data: Dat, weight, config_type
using NBodyIPs: match_dictionary

using Base.Threads

import Base: kron


export get_basis, regression,
       fiterrors, scatter_data,
       print_fiterrors,
       observations, get_lsq_system,
       regularise, table

Base.norm(F::JVecsF) = norm(norm.(F))

# components of the stress (up to symmetry)
const _IS = SVector(1,2,3,5,6,9)

"""
`mutable struct IpLsqSys`: type storing all information to perform a
LSQ fit for an interatomic potential. To assemble the LsqSys use
```
dot(data, basis)
```
"""
mutable struct LsqSys
   data::Vector{Dat}
   basis::Vector{NBodyFunction}
   Iord::Vector{Vector{Int}}
   Ψ::Matrix{Float64}
end


config_types(lsq::LsqSys) = unique(config_type.(lsq.data))

"""
Take a basis and split it into individual body-orders.
"""
function split_basis(basis::AbstractVector{NBodyFunction})
   # get the types of the individual basis elements
   tps = typeof.(basis)
   Iord = Vector{Int}[]
   Bord = Any[]
   for tp in unique(tps)
      # find which elements of basis have type `tp`
      I = find( tp .== tps )
      push!(Iord, I)
      push!(Bord, [b for b in basis[I]])
   end
   return Bord, Iord
end


kron(d::Dat, B::Vector{NBodyFunction}) = dot(d, split_basis(B)...)

# ------- fill the LSQ system, i.e. evaluate basis at data points -------
function kron(d::Dat, Bord::Vector, Iord::Vector{Vector{Int}})
   len = length(d)
   at = Atoms(d)

   # allocate (sub-) matrix of basis functions
   lenY = 1 # energy
   if (forces(d) != nothing) && len > 1
      lenY += 3*len
   end
   if virial(d) != nothing
      lenY += length(_IS)
   end
   Ψ = zeros(lenY, sum(length.(Bord)))

   # energies
   i0 = 0
   for n = 1:length(Bord)
      Es = energy(Bord[n], at)
      Ψ[i0+1, Iord[n]] = Es
   end
   i0 += 1

   # forces
   if (forces(d) != nothing) && len > 1
      for n = 1:length(Bord)
         Fs = forces(Bord[n], at)
         for j = 1:length(Fs)
            fb_vec = mat(Fs[j])[:]
            Ψ[(i0+1):(i0+length(fb_vec)), Iord[n][j]] = fb_vec
         end
      end
      i0 += 3 * len
   end

   # virial components
   if virial(d) != nothing
      for n = 1:length(Bord)
         Ss = virial(Bord[n], at)
         for j = 1:length(Ss)
            Svec = Ss[j][_IS]
            Ψ[(i0+1):(i0+length(_IS)), Iord[n][j]] = Svec
         end
      end
      i0 += length(_IS)
   end

   return Ψ
end



function kron(data::Vector{TD},  basis::Vector{TB}; verbose=true
         ) where {TD <: Dat, TB <: NBodyFunction}
   # sort basis set into body-orders, and possibly different
   # types within the same body-order (e.g. for matching!)
   Bord, Iord = split_basis(basis)

   # generate many matrix blocks, one for each piece of data
   #  ==> this should be switched to pmap, or @parallel
   if nthreads() == 1
      if verbose
         println("Assembly LSQ in serial")
         LSQ = @showprogress(1.0,
                     [ kron(d, Bord, Iord) for d in data ] )
      else
         LSQ = [ kron(d, Bord, Iord) for d in data ]
      end
   else
      # error("parallel LSQ assembly not implemented")
      if verbose
         println("Assemble LSQ with $(nthreads()) threads")
         p = Progress(length(data))
         p_ctr = 0
         p_lock = SpinLock()
      end
      LSQ = Vector{Any}(length(data))
      tic()
      @threads for n = 1:length(data)
         LSQ[n] = kron(data[n], Bord, Iord)
         if verbose
            lock(p_lock)
            p_ctr += 1
            ProgressMeter.update!(p, p_ctr)
            unlock(p_lock)
         end
      end
      verbose && toc()
   end
   # combine the local matrices into a big global matrix
   nrows = sum(size(block, 1) for block in LSQ)
   Ψ = zeros(nrows, length(basis))
   i0 = 0
   for id = 1:length(data)
      Ψi::Matrix{Float64} = LSQ[id]
      i1 = i0 + size(Ψi,1)
      Ψ[(i0+1):i1, :] = Ψi
      i0 = i1
   end

   return LsqSys(data, basis, Iord, Ψ)
end


function observations(d::Dat)
   # ------- fill the data/observations vector -------------------
   Y = Float64[]
   # energy
   push!(Y, energy(d))
   # forces
   if (forces(d) != nothing) && (length(d) > 1)
      f = forces(d)
      f_vec = mat(f)[:]
      append!(Y, f_vec)
   end
   # virial
   if virial(d) != nothing
      S = virial(d)
      append!(Y, S[_IS])
   end
   return Y
end

function observations(data::AbstractVector{Dat})
   Y = Float64[]
   for d in data
      append!(Y, observations(d))
   end
   return Y
end

observations(lsq::LsqSys) = observations(lsq.data)


# ------- Fix a JLD Bug --------------------------------------

import JLD
struct LsqSysSerializer; data; basis; Iord; Ψ; end

JLD.writeas(lsq::LsqSys) = LsqSysSerializer(lsq.data, lsq.basis, lsq.Iord, lsq.Ψ)

function JLD.readas(lsq::LsqSysSerializer)
   basis = lsq.basis
   # make sure all elements of the same basis group have the
   # same dictionary; the problem is that deserialize is called
   # on each NBody individually which will give multiple types that
   # are different for the compiler but describe the same dictionary.
   for I in lsq.Iord, i = 2:length(I)
      basis[I[i]] = match_dictionary(basis[I[i]], basis[I[1]])
   end
   return  LsqSys(lsq.data, basis, lsq.Iord, lsq.Ψ)
end

# -------------------------------------------------------

using NBodyIPs.Data: config_type

function Base.show(io::Base.TTY, lsq::LsqSys)
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


# ------------ Refactored LSQ Fit Code


# _haskey and _getkey are to simulate named tuples

_haskey(t::Tuple, key) = length(find(first.(t) .== key)) > 0

function _getkey(t::Tuple, key)
   i = find(first.(t) .== key)
   return length(i) > 0 ? t[i[1]][2] : nothing
end

function _getkey_val(t::Tuple, key)
   i = find(first.(t) .== key)
   return length(i) > 0 ? t[i[1]][2] : 1.0
end

# TODO: hack - fix it
_getkey_val(::Void, key) = 1.0

default_weights() = (:E => 30.0, :F => 1.0, :V => 0.3)

function analyse_weights(weights::Union{Void, Tuple})
   # default weights
   w_E, w_F, w_V = last.(default_weights())
   if weights != nothing
      _haskey(weights, :E) && (w_E = _getkey(weights, :E))
      _haskey(weights, :F) && (w_F = _getkey(weights, :F))
      _haskey(weights, :V) && (w_V = _getkey(weights, :V))
   end
   return (:E => w_E, :F => w_F, :V => w_V)
end


function analyse_include_exclude(lsq, include, exclude)
   if include != nothing && exclude != nothing
      error("only one of `include`, `exclude` may be different from `nothing`")
   end
   ctypes = config_types(lsq)
   if include != nothing
      if !issubset(include, ctypes)
         error("`include` can only contain config types that are in the dataset")
      end
      # do nothing - just keep `include` as is to return
   elseif exclude != nothing
      if !issubset(exclude, ctypes)
         error("`exclude` can only contain config types that are in the dataset")
      end
      include = setdiff(ctypes, exclude)
   else
      # both are nothing => keep all config_types
      include = ctypes
   end
   return include
end


"""
`get_lsq_system(lsq; kwargs...) -> Ψ, Y, Ibasis`

Assemble the least squares system + rhs. The `kwargs` can be used to
select a subset of the available data or basis set, and to adjust the
weights by config_type. For more complex weight adjustments, one
can directly modify the `lsq.data[n].w` coefficients.

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

* `exclude`, `include`: arrays of strings of config_types to either
include or exclude in the fit (default: all config types are included)

* `order`: integer specifying up to which body-order to include the basis
in the fit. (default: all basis functions are included)
"""
get_lsq_system(lsq; weights=nothing, config_weights=nothing,
                    exclude=nothing, include=nothing, order = Inf,
                    regulariser = nothing) =
   _get_lsq_system(lsq, analyse_weights(weights), config_weights,
                   analyse_include_exclude(lsq, include, exclude), order,
                   regulariser)

# function barrier for get_lsq_system
function _get_lsq_system(lsq, weights, config_weights, include, order, regulariser)

   Y = observations(lsq)
   W = zeros(length(Y))
   # energy, force, virial weights
   w_E, w_F, w_V = last.(weights)
   # reference energy => we assume the first basis function is 1B
   E0 = lsq.basis[1]()

   # assemble the weight vector
   idx = 0
   for d in lsq.data
      len = length(d)
      # weighting factor due to config_type
      w_cfg = _getkey_val(config_weights, config_type(d))
      # weighting factor from dataset
      w = weight(d) * w_cfg

      # keep going through all data, but set the weight to zero if this
      # one is to be excluded
      if !(config_type(d) in include)
         w = 0.0
      end

      # energy
      W[idx+1] = w * w_E
      idx += 1
      # and while we're at it, subtract E0 from Y
      Y[idx] -= E0 * len

      # forces
      if (forces(d) != nothing) && (len > 1)
         W[(idx+1):(idx+3*len)] = w * w_F
         idx += 3*len
      end
      # virial
      if virial(d) != nothing
         W[(idx+1):(idx+length(_IS))] = w * w_S
         idx += length(_IS)
      end
   end
   # double-check we haven't made a mess :)
   @assert idx == length(W) == length(Y)

   # find the zeros and remove them => list of data points
   Idata = find(W .!= 0.0) |> sort

   # find all basis functions with the required body-order
   # (note we also remove the B1, which is assumed to be at index 1)
   Ibasis = find(1 .< bodyorder.(lsq.basis) .<= order) |> sort

   # take the appropriate slices of the data and basis
   Y = Y[Idata]
   W = W[Idata]
   Ψ = lsq.Ψ[Idata, Ibasis]

   # now rescale Y and Ψ according to W => Y_W, Ψ_W; then the two systems
   #   \| Y_W - Ψ_W c \| -> min  and (Y - Ψ*c)^T W (Y - Ψ*x) => MIN
   # are equivalent
   W .= sqrt.(W)
   Y .*= W
   scale!(W, Ψ)

   if any(isnan, Ψ) || any(isnan, Y)
      error("discovered NaNs - something went horribly wrong!")
   end

   # regularise
   if regulariser != nothing
      P = regulariser(lsq.basis[Ibasis])
      Ψ, Y = regularise(Ψ, Y, P)
   end

   # this should be it ...
   return Ψ, Y, Ibasis
end


function regularise(Ψ, Y, P::Matrix)
   @assert size(Ψ,2) == size(P,2)
   return vcat(Ψ, P), vcat(Y, zeros(size(P,1)))
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
function regression(lsq;
                    verbose = true,
                    kwargs...)
   # TODO
   #  * regulariser
   #  * stabstyle or stabiliser => think about this carefully!

   # apply all the weights, get rid of anything that isn't needed or wanted
   # in particular subtract E0 from the energies and remove B1 from the
   # basis set
   Y, Ψ, Ibasis = get_lsq_system(lsq; kwargs...)

   # QR factorisation
   verbose && println("solve $(size(Ψ)) LSQ system using QR factorisation")
   QR = qrfact(Ψ)
   verbose && @show cond(QR[:R])
   # back-substitution to get the coefficients # same as QR \ ???
   c = QR \ Y

   # check error on training set: i.e. the naive errors using the
   # weights that are provided
   if verbose
      z = Ψ * c - Y
      rel_rms = norm(z) / norm(Y)
      verbose && println("naive relative rms error on training set: ", rel_rms)
   end

   # now add the 1-body term and convert this into an IP
   # (I assume that the first basis function is the 1B function)
   @assert bodyorder(lsq.basis[1]) == 1
   return NBodyIP(lsq.basis, [1.0; c])
end


function NBodyIP(lsq::LsqSys, c::Vector, Ibasis::Vector{Int})
   if !(1 in Ibasis)
      @assert bodyorder(lsq.basis[1]) == 1
      c = [1.0; c]
      Ibasis = [1; Ibasis]
   end
   return NBodyIP(lsq.basis[Ibasis], c)
end


# =============== Analysis Fitting Errors ===========

_fiterrsdict() = Dict("E-RMS" => 0.0, "F-RMS" => 0.0, "E-MAE" => 0.0, "F-MAE" => 0.0)

struct FitErrors
   errs::Dict{String, Dict{String,Float64}}
   nrms::Dict{String, Dict{String,Float64}}
end

function fiterrors(lsq, c, Ibasis; include=nothing, exclude=nothing)
   include = analyse_include_exclude(lsq, include, exclude)

   E0 = lsq.basis[1]()
   # create the dict for the fit errors
   errs = Dict( "set" => _fiterrsdict() )
   # count number of configurations
   num = Dict("set" => Dict("E" => 0, "F" => 0))
   obs = Dict("set" => _fiterrsdict())

   for ct in include
      errs[ct] = _fiterrsdict()
      num[ct] = Dict("E" => 0, "F" => 0)
      obs[ct] = _fiterrsdict()
   end

   idx = 0
   for d in lsq.data
      ct = config_type(d)
      E_data, F_data, len = energy(d), forces(d), length(d)
      E_data -= E0 * len

      if !(ct in include)
         idx += 1 # E
         if F_data != nothing
            idx += 3 * len
         end
         if virial(d) != nothing
            idx += length(_IS)
         end
         continue
      end


      # ----- Energy error -----
      E_fit = dot(lsq.Ψ[idx+1, Ibasis], c)
      Erms = (E_data - E_fit)^2 / len^2
      Emae = abs(E_data - E_fit) / len
      errs[ct]["E-RMS"] += Erms
      obs[ct]["E-RMS"] += (E_data / len)^2
      errs[ct]["E-MAE"] += Emae
      obs[ct]["E-MAE"] += abs(E_data / len)
      num[ct]["E"] += 1
      errs["set"]["E-RMS"] += Erms
      obs["set"]["E-RMS"] += (E_data / len)^2
      errs["set"]["E-MAE"] += Emae
      obs["set"]["E-MAE"] += abs(E_data / len)
      num["set"]["E"] += 1
      idx += 1

      # ----- Force error -------
      if F_data != nothing
         f_data = mat(F_data)[:]
         f_fit = lsq.Ψ[(idx+1):(idx+3*len), Ibasis] * c
         Frms = norm(f_data - f_fit)^2
         Fmae = norm(f_data - f_fit, 1)
         errs[ct]["F-RMS"] += Frms
         obs[ct]["F-RMS"] += norm(f_data)^2
         errs[ct]["F-MAE"] += Fmae
         obs[ct]["F-MAE"] += norm(f_data, 1)
         num[ct]["F"] += 3 * len
         errs["set"]["F-RMS"] += Frms
         obs["set"]["F-RMS"] += norm(f_data)^2
         errs["set"]["F-MAE"] += Fmae
         obs["set"]["F-MAE"] += norm(f_data, 1)
         num["set"]["F"] += 3 * len
         idx += 3 * len
      end

      # skip the stresses
      if virial(d) != nothing
         idx += length(_IS)
      end
   end

   # NORMALISE
   for key in keys(errs)
      nE = num[key]["E"]
      nF = num[key]["F"]
      errs[key]["E-RMS"] = sqrt(errs[key]["E-RMS"] / nE)
      obs[key]["E-RMS"] = sqrt(obs[key]["E-RMS"] / nE)
      errs[key]["F-RMS"] = sqrt(errs[key]["F-RMS"] / nF)
      obs[key]["F-RMS"] = sqrt(obs[key]["F-RMS"] / nF)
      errs[key]["E-MAE"] = errs[key]["E-MAE"] / nE
      obs[key]["E-MAE"] = obs[key]["E-MAE"] / nE
      errs[key]["F-MAE"] = errs[key]["F-MAE"] / nF
      obs[key]["F-MAE"] = obs[key]["F-MAE"] / nF
   end

   return FitErrors(errs, obs)
end


function table(errs::FitErrors; relative=false)
   if !relative
      println(errs::FitErrors)
      return
   end
   print("-------------------------------------------------\n")
   print("             ||       RMSE     ||       MAE      \n")
   print(" config type ||  E [%] | F [%] ||  E [%] | F [%] \n")
   print("-------------||--------|-------||--------|-------\n")
   s_set = ""
   nrms = errs.nrms
   for (key, D) in errs.errs
      nrm = nrms[key]
      lkey = min(length(key), 11)
      s = @sprintf(" %11s || %1.4f | %1.3f || %1.4f | %1.3f \n",
               key[1:lkey], D["E-RMS"]/nrm["E-RMS"], D["F-RMS"]/nrm["F-RMS"],
                            D["E-MAE"]/nrm["E-MAE"], D["F-MAE"]/nrm["F-MAE"])
      if key == "set"
         s_set = s
      else
         print(s)
      end
   end
   print("-------------||--------|-------||--------|-------\n")
   print(s_set)
   print("-------------------------------------------------\n")
end


function Base.show(io::Base.TTY, errs::FitErrors)
   print(io, "-------------------------------------------------------\n")
   print(io, "             ||       RMSE        ||        MAE       \n")
   print(io, " config type ||  E [eV] | F[eV/A] ||  E [eV] | F[eV/A] \n")
   print(io, "-------------||---------|---------||---------|---------\n")
   s_set = ""
   for (key, D) in errs.errs
      lkey = min(length(key), 11)
      s = @sprintf(" %11s || %1.5f | %1.5f || %1.5f | %1.5f \n",
               key[1:lkey], D["E-RMS"], D["F-RMS"], D["E-RMS"], D["E-MAE"])
      if key == "set"
         s_set = s
      else
         print(io, s)
      end
   end
   print(io, "-------------||---------|---------||---------|---------\n")
   print(io, s_set)
   print(io, "-------------------------------------------------------\n")
end



#
# """
# computes the maximum force over all configurations
# in `data`
# """
# function max_force(b, data)
#    out = 0.0
#    for d in data
#       f = forces(b, Atoms(d))
#       out = max(out, maximum(norm.(f)))
#    end
#    return out
# end
#
# function normalize_basis!(B, data)
#    for b in B
#       @assert (length(b.c) == 1)
#       maxfrc = max_force(b, data)
#       if maxfrc > 0
#          b.c[1] /= maxfrc
#       end
#       if 0 < maxfrc < 1e-8
#          warn("encountered a very small maxfrc = $maxfrc")
#       end
#    end
#    return B
# end
#


# function scatter_data(IP, data)
#    E_data = Float64[]
#    E_fit = Float64[]
#    F_data = Float64[]
#    F_fit = Float64[]
#    for d in data
#       at = Atoms(d)
#       len = length(at)
#       push!(E_data, energy(d) / len)
#       append!(F_data, mat(forces(d))[:])
#       push!(E_fit, energy(IP, at) / len)
#       append!(F_fit, mat(forces(IP, at))[:])
#    end
#    return E_data, E_fit, F_data, F_fit
# end
