
"""
`module Lsq`

This sub-modules contains code to convert a database `LsqDB` into
linear LSQ systems, apply suitable weights, solve the resulting systems
and return an IP (e.g. an NBodyIP).

The ordering of observations (entries of Y, row-indices of Psi) is given by
the loop ordering
```
for ct in configtypes
   for dt in datatypes  # observation types
      for dat in db.data_groups[ct]
      # or: db.kron_groups[ct][dt][:,:,...]
```
E.g., if the configtypes are "E", "F", then
```
Y = [E, E, E, E...., F, F, F, F, ...]
```
"""
module Lsq

import JuLIP, IPFitting, StatsBase

using StaticArrays
using JuLIP: vecs
using JuLIP: AbstractCalculator, Atoms
using JuLIP.Potentials: OneBody
using JuLIP.MLIPs: SumIP
using IPFitting: Dat, LsqDB, weighthook, observations,
                      observation, hasobservation, eval_obs, vec_obs
using IPFitting.Data: configtype
using IPFitting.DB: dbpath, _nconfigs, matrows
using LinearAlgebra: lmul!, Diagonal, qr, qr!, cond, norm, svd, I, pinv, mul!, cholesky
using InteractiveUtils: versioninfo
using LowRankApprox
using JuLIP.Utils
using JuLIP: JVecF
using PyCall
using IterativeSolvers


const Err = IPFitting.Errors

export lsqfit, onb, stepwise_fit

"""
`collect_observations(db::LsqDB, weights::Dict, Vref) -> Y, W, Icfg

Construct vectors of observations and weights from the lsq system specified
by `db`. The observations do not have the weights applied yet.

Weights are specified as follows:
```
   wh = weighthook(obskey, cfg)
   w = w1 * wh
```
where `w1` is obtained from the `weights` dictionary.
```
weights = Dict( "default" => Dict( "E" => 30.0, "F" => 1.0, "V" => 0.3 ),
                "bulk" => Dict( "E" => 100.0, "F" = 20.0, "V" => 0.0 ),
                ... )
```
To compute the weights, the code first looks for the `configtype` entry. If
it doesn't exist, then it reverts to the `default` entry. The weights obtained
from the `weights` dictionary is then multiplied by a

- The above procedure can be overwritten by an entry in `cfg.info["W"]`;
in this case the re-weighting via the `weighthook` is ignored as well.
- If `configtype(cfg) ∈ weights["ignore"]` then this configuration is skipped.
- If `obskey ∈ weights["ignore"]` then this observation is skipped
"""
function collect_observations(db::LsqDB,
                              weights::Dict,
                              Vref)#, scal_wgs = false )

   nrows = size(db.Ψ, 1)  # total number of observations if we collect
                          # everything in the database
   Y = zeros(Float64, nrows)
   W = zeros(Float64, nrows)
   Icfg = zeros(Int, nrows)
   for (obskey, dat, icfg) in observations(db)  # obskey ∈ {"E","F",...}; d::Dat
      if configtype(dat) in weights["ignore"] || obskey in weights["ignore"]
         continue
      end

      irows = matrows(dat, obskey)
      ct = configtype(dat)

      # ----- Observation ------
      obs = observation(obskey, dat)
      # If Vref != nothing the it is a calculator and we can subtract
      # its value => this allows us to fit from a reference potential
      if haskey(dat.info, "Vref")
         obs = obs - dat.info["Vref"][obskey]
      elseif Vref != nothing
         obs = obs - vec_obs(obskey, eval_obs(obskey, Vref, dat))
      end
      Y[irows] .= obs
      Icfg[irows] .= icfg
      # ------ Weights --------
      # modify the weights from extra information in the dat structure
      W[irows] .= _get_weights(weights,
                               weighthook(obskey, dat),
                               dat, obskey, obs)#, scal_wgs)
   end
   return Y, W, Icfg
end

#f_w(fi, fm; A=0.5, B=0.5, f0=0.5) = (A + (B * f0 * log(1 + fi/f0 + fm/f0)))^(-1.0)
f_w(fi, fm) = (fi + fm)^(-0.5)


"""
see documentation in `collect_observations`
"""
function _get_weights(weights, wh, dat, obskey, o)
   # initialise the weight to 0.0. If this isn't overwritten then it means
   # we will ignore this observation.
   w = 0.0

   # if there is a "W" entry in dat.info then this means all defaults
   # are over-written and we use this weight directly as is
   if haskey(dat.info, "W")
      # now check that this observation type (E, F, V) exists in dat.D["W"]
      # if yes return the corresponding weight, if not we ignore that
      # observation.
      if haskey(dat.info["W"], obskey)
         w = dat.info["W"][obskey]
      end
   else
      # if no "W" dict exists extract the correct weights from
      # the `weights` dictionary
      # first check what the current configtype has a separate entry
      # and otherwise use the default
      ct = configtype(dat)
      if haskey(weights, ct)
         cfgkey = ct
      else
         cfgkey = "default"
      end
      # now check whether there is an entry for the current observation  type
      if haskey(weights[cfgkey], obskey)
         w = weights[cfgkey][obskey] * wh
      end
   end

   if obskey in ["E", "V"]
   # transform the weights into a vector (if necessary) and return
      if length(w) == 1
         return w * ones(length(o))
      elseif length(w) == length(o)
         return w
      end
   elseif obskey == "F"
      wf = w * ones(length(o))#* f_w.(abs.(dat.D["F"]), mean(abs.(dat.D["F"]))) 
      return wf
   end
   error("_get_weights: length(w) is neither 1 nor length(o)?!?!?")
end

function _regularise!(Ψ::Matrix{T}, Y::Vector{T}, basis, regularisers;
                      verbose=false, Ibasis = :) where {T}
   # TODO: incorporate Ibasis again !
   # assemble the regularisers
   Ψreg = Matrix{Float64}(undef, 0, size(Ψ, 2))
   Yreg = Float64[]
   for reg in regularisers
      if reg isa Matrix
         P = reg
         q = zeros(size(P, 1))
      else
         P, q = Matrix(reg, basis; verbose=verbose)
      end
      Ψreg = vcat(Ψreg, P)
      Yreg = vcat(Yreg, q)
   end
   # check they all have the correct size
   return vcat(Ψ, Ψreg), vcat(Y, Yreg)
end
#using DelimitedFiles

function _forceprecon(db, Y, weights)
   if !haskey(weights, "precon")
      return db.Ψ, Y
   end
   if isempty(weights["precon"])
      return db.Ψ, Y
   end

   Ψ = copy(db.Ψ)

   for (obskey, dat, icfg) in observations(db)  # obskey ∈ {"E","F",...}; d::Dat
      if obskey != "F"; continue; end
      ct = configtype(dat)
      if !haskey(weights["precon"], ct) && !haskey(weights["precon"], "default"); continue; end  ###### TOBE changed if we want to precondition only certain config_types

      #preconfun = ["precon"][ct]
      if haskey(weights["precon"], ct)
         P = weights["precon"][ct](dat.at)
      elseif haskey(weights["precon"], "default")
         P = weights["precon"]["default"](dat.at)
      end
      #   E(R + U) - E(R) - DE(R).U ~ U' D^2E(R) U ~ U' P U -> sqrt(U' P U)
      # DE(R+U) - DE(R) ~ D^2E(R) U
      # instead : DE(R+U) - DE(R) ~ P U
      #         || P^{-1/2} (DE(R+U) - DE(R)) || ~ || P^{1/2} U || = sqrt(U' P U)
      #         || P^{-1} (DE(R+U) - DE(R)) || ~ || U ||  => PROBABLY NOT THIS!!!!
      # all we need is
      #     sqrt(P)^2 = P
      #    || sqrt(P)^{-1} F ||^2 = F' sqrt{P}^{-2} F = F' P^{-1} F
      # Option 1: sqrt(P) = Q sqrt(Λ) Q'
      # Option 2: sqrt(P) = chol(P).L
      if !haskey(weights["precon"],"_solver")
         solver = "chol"
      else
         solver = weights["precon"]["_solver"]
      end
      if solver == "chol"
         rtP = cholesky(P, Val(false)).L
      elseif solver == "sqrt"
         rtP = sqrt(P)
      else
         error("unknown root solver")
      end

      irows = matrows(dat, obskey)
      Y[irows] = rtP \ Y[irows]
      Ψ[irows, :] = rtP \ Ψ[irows, :]
   end
   @info("USING FORCE PRECONDITIONER")
   return Ψ, Y
end

"""
`get_lsq_system(db::LsqDB; kwargs...) -> Ψ, Y`

Assemble the least squares system + rhs (observations). The `kwargs` can be used
to select a subset of the available data or basis set, and to adjust the
weights by config_type and observation type. See `lsqfit` for a list
of keyword arguments. Allowed/required kwargs are `verbose, weights,
E0, Vref, Ibasis, Itrain, regularisers`.
"""
@noinline function get_lsq_system(db::LsqDB;
                                  verbose=true,
                                  Ibasis = :,
                                  Itrain = :,
                                  E0 = nothing,
                                  Vref = OneBody(E0),
                                  weights = nothing,
                                  #scal_wgs = false,
                                  regularisers = [])
   weights = _fix_weights!(weights)

   # get the observations vector and the weights vector
   # the Vref potential is subtracted from the observations
   Y, W, Icfg = collect_observations(db, weights, Vref)#, scal_wgs)

   # check for NaNs
   any(isnan, Y) && error("NaN detected in observations vector")
   any(isnan, W) && error("NaN detected in weights vector")

   # NON-DIAGONAL ADDITIONAL WEIGHTS
   #  If the user supplies a force preconditioner, then we premultiply
   #  Ψ and Y with another weight function that modifies everything
   Ψ, Y = _forceprecon(db, Y, weights)

   # convert : into a vector or list
   if Itrain == : ; Itrain = 1:length(Y); end

   # get the right slice of the big fat huge enourmous LSQ matrix
   # the columns are just the basis functions
   Icols = Ibasis
   # the rows are those which have (a) a non-zero weight
   # and (b) the configuration index is part of the training set
   Irows = intersect( findall(W .!= 0) |> sort,
                      findall(in.(Icfg, Ref(Itrain))) )

   # we make this a view but make sure to copy it before applying the weights
   Ψ = @view Ψ[Irows, Icols]
   Y = Y[Irows]
   W = W[Irows]

   if any(isnan, Ψ)
      @error("discovered NaNs in LSQ system matrix")
      naninds = findall(isnan.(Ψ))
      n = min(length(naninds), 10)
      @show length(naninds)
      @show naninds[1:n]
      error("discovered NaNs in LSQ system matrix")
   end

   # regularise
   if !isempty(regularisers)
      Ψ, Y = _regularise!(Ψ, Y, db.basis, regularisers; verbose=verbose, Ibasis=Ibasis)
      # adjust the weight-vector to the new system size
      append!(W, ones(length(Y)-length(W)))
   else
      Ψ = collect(Ψ)
   end

   # now rescale Y and Ψ according to W => Y_W, Ψ_W; then the two systems
   #   \| Y_W - Ψ_W c \| -> min  and \| W (Y - Ψ*c)^T \| -> min
   # are equivalent
   Y = Diagonal(W) * Y
   lmul!(Diagonal(W), Ψ)

   # this should be it ...
   return Ψ, Y
end

_show_free_mem() =
      @info("Free Memory: ≈ $(round(Sys.free_memory()*1e-9, digits=2)) GB")


"""
`lsqfit(db; kwargs...) -> IP, errs`

Given the pre-computed least-squares system `db` (cf `LsqDB`) setup a least
squares system, compute the solution and return an interatomic potential
`IP` as well as the associated error estimates `errs`.

## Keyword Arguments:

* `E0` (required unless `Vref` is provided) : energy of the `OneBody` term
* `Vref` : a reference potential from which to start the fit.
* `weights` (required) : A dictionary specifying the weights for different
types of configurations, see `?collect_observations` for how this is specified.
* `Ibasis` : indices of basis functions to be used in the fit, default is `:`
* `verbose` : true or false
* `solver` : -experimental, still need to  write the docs for this-
* `regularisers` : a list of regularisers to be added to the lsq functional.
Each regulariser `R` can be of an arbtirary type but this type must implement the
conversion to matrix
```
Matrix(R, basis; verbose={true,false})
```
If `Areg = Matrix(R, basis)` then this corresponds to adding
`|| Areg * x ||²` to the least squares functional.
* `deldb = false` : experimental - if `true` then the lsq matrix in the `db` is
deleted after assembling the weighted lsq system
* `error_table = false` : experimental - if `true` then the lsq matrix is used to
assemble the fitting errors.
* `saveqr = nothing` : experimental - if `saveqr` is a Dict then the factors
of the QR factorisation are stored in it for future use outside of this
function

## Return types

* `IP`: whatever `JuLIP.MLIPs.combine` returns
* `fitinfo::Dict`: adds all sort of information about the fit
"""

function do_qr!(db, Ψ, Y)
   @info("Performing QR decomposition")
   GC.gc();
   qrΨ = qr!(Ψ)

   saveqr = Dict()
   saveqr["Q"] = Matrix(qrΨ.Q)
   saveqr["R"] = Matrix(qrΨ.R)
   saveqr["Y"] = Y

   db = LsqDB(db.basis, db.configs, db.Ψ, db.dbpath, saveqr)

   flush_qr(db)

   return db
end

"""
`stepwise_fit()`

Function fitting first dB1, than using the fit as a baseline fitting dB2
Can be used for example to fit a 2B first, and than an ACE afterwards, not fitting them together

Returns sumIP of Vref, basis1 and basis2. 
"""

@noinline function stepwise_fit(dB1::LsqDB, 
                              dB2::LsqDB;
                              solver1 = Dict("solver" => :qr),    
                              solver2 = solver1, 
                              Vref = nothing,
                              weights = nothing,
                              error_table = false,
                              verbose = true)
   IP_1, lsqinfo1 = lsqfit(dB1, 
                           solver=solver1,
                           Vref = Vref,
                           weights = weights,
                           error_table = error_table,
                           verbose=verbose)
   IP_2, lsqinfo2 = lsqfit(dB2, 
                           solver=solver2,
                           Vref = IP_1,
                           weights = weights,
                           error_table = error_table,
                           verbose=verbose)
   return IP_2, lsqinfo2
end

@noinline function lsqfit(db::LsqDB;
                solver=Dict("solver" => :qr), verbose=true,
                Ibasis = :,
                Itrain = :,
                #Itest = nothing,
                E0 = 0.0,
                Vref = nothing, # OneBody(E0),
                weights = nothing,
                sigmas = nothing,
                regularisers = [],
                deldb = false,
                error_table = false,
                #lasso = false,
                saveqr = nothing,
                #scal_wgs = false,
                kwargs...)

   if !isnothing(weights) && !isnothing(sigmas)
      @error("cannot specify both weights and sigmas")
   elseif !isnothing(sigmas)
      weights = sigmas
      for (config, config_weights) in weights
         for (key, weight) in config_weights
            config_weights[key] = 1.0/config_weights[key]
         end
      end
   end
   weights = _fix_weights!(weights)
   Jbasis = ((Ibasis == Colon()) ? (1:length(db.basis)) : Ibasis)

   verbose && @info("assemble lsq system")
   verbose && _show_free_mem()
   Ψ, Y = get_lsq_system(db; verbose=verbose, Vref=Vref,
                             Ibasis=Ibasis, Itrain = Itrain,
                             weights = weights,
                             regularisers = regularisers)
                             #scal_wgs = scal_wgs)
   verbose && _show_free_mem()
   if deldb
      @info("Deleting database - `db` can no longer be saved to disk")
      db.Ψ = Matrix{Float64}(undef, 0,0)
      db.dbpath = ""
   end
   GC.gc()
   verbose && _show_free_mem()

   if !haskey(solver, "P")
      P = Diagonal(I, length(Ψ[1,:]))
   else
      @info("Using Laplacian Preconditioning (P)")
      P = solver["P"]
   end

   D_inv = pinv(P)
   mul!(Ψ,Ψ,D_inv)
   
   if haskey(solver, "P")
      delete!(solver, "P")
   end

   #κ, p_1, int_order = 0.0, 0.0, 0.0
   @info("Size of least squares system: $(size(Ψ))")
   if solver["solver"] == :qr
      @info("Using QR")
      qrΨ = qr!(Ψ)
      κ = cond(qrΨ.R)
      GC.gc();
      verbose && @info("cond(R) = $(cond(qrΨ.R))")
      c = qrΨ \ Y
      rel_rms = norm(qrΨ.Q * (qrΨ.R * c) - Y) / norm(Y)

      if saveqr isa Dict
         saveqr["Q"] = qrΨ.Q
         saveqr["R"] = qrΨ.R
         saveqr["Y"] = Y
      end

      qrΨ = nothing
   elseif solver["solver"] == :svd
      @assert haskey(solver, "svd_ndiscard")
      svd_ndiscard = solver["svd_ndiscard"]
      @info("Using SVD: svd_ndiscard=$(svd_ndiscard)")
      F = svd(Ψ)
      c = F.V[:,1:(end-svd_ndiscard)] * (Diagonal(F.S[1:(end-svd_ndiscard)]) \ (F.U' * Y)[1:(end-svd_ndiscard)])
      rel_rms = norm(Ψ * c - Y) / norm(Y)
   elseif solver["solver"] == :rrqr
      @assert haskey(solver, "rrqr_tol")
      rrqr_tol = solver["rrqr_tol"]
      @info("Using RRQR: rrqr_tol=$(rrqr_tol)")
      qrΨ = pqrfact(Ψ, rtol=rrqr_tol)
      verbose && @info("cond(R) = $(cond(qrΨ.R))")
      c = qrΨ \ Y
      rel_rms = norm(Ψ * c - Y) / norm(Y)
   elseif solver["solver"] == :lsqr
      @assert haskey(solver, "lsqr_damp")
      @assert haskey(solver, "lsqr_atol")
      lsqr_damp = solver["lsqr_damp"]
      lsqr_atol = solver["lsqr_atol"]
      lsqr_conlim = haskey(solver,"lsqr_conlim")   ? solver["lsqr_conlim"]  : 1e8
      lsqr_maxiter = haskey(solver,"lsqr_maxiter") ? solver["lsqr_maxiter"] : maximum(size(Ψ))
      lsqr_verbose = haskey(solver,"lsqr_verbose") ? solver["lsqr_verbose"] : false
      if haskey(solver, "c_init")
         @info("Using c_init")
         c_init = solver["c_init"]
      else
         c_init = zeros(length(Ψ[1,:]))
      end
      @info("Using LSQR: damp=$(lsqr_damp), atol=$(lsqr_atol), conlim=$(lsqr_conlim), maxiter=$(lsqr_maxiter)")
      c, ch = lsqr!(c_init, Ψ, Y, damp=lsqr_damp, atol=lsqr_atol, conlim=lsqr_conlim,
                          maxiter=lsqr_maxiter, verbose=lsqr_verbose, log=true)

      if ch[:anorm][end] > lsqr_atol
         print("\n\n")
         @warn "LSQR did not converge (atol=$(lsqr_atol) but anorm=$(ch[:anorm][end])). Consider increasing conlim or maxiter."
         print("\n\n")
      end

      rel_rms = norm(Ψ * c - Y) / norm(Y)
   elseif solver["solver"] == :brr
      BRR = pyimport("sklearn.linear_model")["BayesianRidge"]
      @assert haskey(solver, "brr_tol")
      if !haskey(solver, "brr_fit_intercept")
         fit_intercept=true
      else
         fit_intercept=solver["brr_fit_intercept"]
      end
      brr_tol = solver["brr_tol"]
      @info("Using BRR: brr_tol=$(brr_tol)")

      clf = BRR(tol=brr_tol, fit_intercept=fit_intercept, normalize=true, compute_score=true)
      clf.fit(Ψ, Y)

      c = clf.coef_
      alpha = clf.alpha_
      lambda = clf.lambda_
      score = clf.scores_[end]
      @info("alpha=$(alpha), lambda=$(lambda), score=$(score)")

      rel_rms = norm(Ψ * c - Y) / norm(Y)
   elseif solver["solver"] == :ard
      ARD = pyimport("sklearn.linear_model")["ARDRegression"]
      @assert haskey(solver, "ard_threshold_lambda")
      @assert haskey(solver, "ard_tol")
      if !haskey(solver, "ard_fit_intercept")
         fit_intercept=true
      else
         fit_intercept=solver["ard_fit_intercept"]
      end
      ard_threshold_lambda = solver["ard_threshold_lambda"]
      ard_tol = solver["ard_tol"]
      @info("Using ARD: ard_tol=$(ard_tol), ard_threshold_lambda=$(ard_threshold_lambda)")

      clf = ARD(threshold_lambda = ard_threshold_lambda, tol=ard_tol, fit_intercept=fit_intercept, normalize=true, compute_score=true)
      clf.fit(Ψ, Y)

      c = clf.coef_
      alpha = clf.alpha_
      lambda = clf.lambda_
      score = clf.scores_[end]

      non_zero_ind = findall(x -> x != 0.0, c)

      @info("Fit complete: keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(c), digits=2)*100)%)")
      @info("score=$(score)")

      rel_rms = norm(Ψ * c - Y) / norm(Y)
   elseif solver["solver"] == :xgboost
      XGBoost = pyimport("xgboost")["XGBRegressor"]

      model = XGBoost(booster="gblinear")
      model.fit(Ψ, Y)

      c = model.coef_

      rel_rms = norm(Ψ * c - Y) / norm(Y)
   else
      error("unknown `solver` in `lsqfit`")
   end

   c = D_inv * c

   Ψ = nothing
   GC.gc()
   verbose && _show_free_mem()

   if verbose
      @info("Relative RMSE on training set: $rel_rms")
   end

   IP = JuLIP.MLIPs.combine(db.basis, c)
   if Vref != nothing
      IP = SumIP(Vref, IP)
   end

   infodict = asm_fitinfo(db, IP, c, Ibasis, weights,
                          Vref, E0, regularisers, verbose,
                          #Itrain, Itest, 
                          error_table)
   GC.gc()
   return IP, merge(infodict, solver)
end


function asm_fitinfo(db, IP, c, Ibasis, weights,
                     Vref, E0, regularisers, verbose,
                     #Itrain = :, Itest = nothing, 
                     error_table=true)
   if Ibasis isa Colon
      Jbasis = collect(1:length(db.basis))
   else
      Jbasis = Ibasis
   end

   cfgtypes = setdiff(unique(configtype.(db.configs)), weights["ignore"])

   # compute errors TODO: still need to fix this!
   # if error_table
   #    verbose && @info("Assemble errors table")
   #    @warn("new error implementation... redo this part please ")
   #    errs = Err.lsqerrors(db, c, Jbasis;
   #             cfgtypes=cfgtypes, Vref=Vref, Icfg=Itrain)
   #    if Itest != nothing
   #       errtest = Err.lsqerrors(db, c, Jbasis;
   #                cfgtypes=cfgtypes, Vref=Vref, Icfg=Itest)
   #    else
   #       errtest = Dict()
   #    end
   # end
   if error_table
      @info("Assembling Error Table")
      Err.add_fits_serial!(IP, db.configs, fitkey="IP")
      rmse_, rmserel_ = Err.rmse(db.configs; fitkey="IP");
      errs = Dict("rmse" => rmse_,
                  "relrmse" => rmserel_)
   end
   # --------------------------------------------------------------------
   # ASSEMBLE INFO DICT
   # --------------------------------------------------------------------
   verbose && @info("Assemble Information about the fit")

   # Julia Version Info
   iob = IOBuffer()
   versioninfo(iob)
   #juliainfo = String(take!(iob))

   infodict = Dict(#"solver" => String(solver[1]),
                   "E0"     => E0,
                   "Ibasis" => Vector{Int}(Jbasis),
                   "c"      => c,
                   "dbpath" => dbpath(db),
                   "weights" => weights,
                   "regularisers"  => Dict.(regularisers),
                   #"juliaversion"  => juliainfo,
                   "errors" => Dict(),
                   "IPFitting_version" => "na", # get_pkg_info("IPFitting"),
                  )
   # TODO: fix IPFitting_version retrieval

   if error_table
      infodict["errors"] = errs
      #infodict["errtest"] = errtest
   end
   # --------------------------------------------------------------------
   return infodict
end

import Pkg

# function get_pkg_info(pkg::AbstractString)
#     pkgs = [Pkg.API.check_package_name(pkg)]
#     ctx = Pkg.Types.Context()
#     Pkg.API.project_resolve!(ctx.env, pkgs)
#     Pkg.API.project_deps_resolve!(ctx.env, pkgs)
#     Pkg.API.manifest_resolve!(ctx.env, pkgs)
#     Pkg.API.ensure_resolved(ctx.env, pkgs)
#     i = Pkg.Display.status(ctx, pkgs)[1]
#     return Dict("name" => i.name,
#                 "uuid" => string(i.uuid),
#                 "ver" => string(i.new.ver))
# end


_fix_weights!(::Nothing) = _fix_weights!(Dict{String, Any}())

function _fix_weights!(weights::Dict)
   if !haskey(weights, "ignore")
      weights = Dict{Any,Any}(weights)
      weights["ignore"] = String[]
   end
   if !haskey(weights, "default")
      weights["default"] = Dict("E" => 1.0, "F" => 1.0, "V" => 1.0)
   end
   return weights
end



"""
`function reglsq(; kwargs...)`

Solves the problem
```
   min || Γ * x || subj to || A * x - z || ≦ tol
```

### Arguments:
* `Γ, R, z, tol` : specify the problem to be solved
* Instead of `tol` one may pass the kwargs `η0, τ`, which will set `tol = sqrt(τ^2-η0^2)`
* `abstol = 1e-2, reltol = 1e-2` : termination tolerance
* `λinit = 1e-3, λmax = 1e2` : search parameters
"""
function reglsq( ;
                Γ = nothing,
                R = nothing,
                y = nothing,
                η0 = nothing,
                τ = nothing,
                tol = sqrt(τ^2 - η0^2),
                abstol = 1e-5,
                reltol = 1e-2,
                λinit = 1e-3,
                λmax = 1e2,
                verbose = true)
   @assert size(Γ, 2) == size(R, 2)
   @assert length(y) == size(R, 1)
   Qc, Rc = qr(collect(Γ))
   Y = [y; zeros(size(Γ, 1))]

   verbose && @info("`reglsq` : solve regularised least squares")
   x = nothing

   function _solve(λ)
      x = qr([ R; λ * Rc ]) \ Y
      return x, norm(R * x - y)
   end

   # Bracketing
   λ1 = 0.0
   λ2 = λinit
   while true
      x, η = _solve(λ2)
      if η > tol
         break
      end
      if λ2 > λmax
         @error("""reglsq reached λmax = $λmax in the bracketing stage.
                  The returned coeffcients are in the interior of the constraint
                  set ||Rx-z|| ≦ tol. If stronger regularisation is required,
                  then please increase the λmax parameter.""")
         return x
      end

      λ1 = λ2
      λ2 *= 10
   end

   verbose && @info("found bracket, starting bisection")

   # Bisection
   while true
      λ = 0.5 * (λ1 + λ2)
      x, η = _solve(λ)
      if tol / (1+reltol) <= η <= tol * (1+reltol)
         # success
         break
      end
      if η > tol
         λ2 = λ
      else
         λ1 = λ
      end
      if abs(λ1 - λ2) < abstol
         @error("""bracket has become too small; returning current solutions,
                   if it is unsatisfactory, then please reduce the
                   `abstol` parameter.""", λ1, λ2, abstol)
         return x
      end
   end

   verbose && @info("found a solution")
   return x
end


function reglsq(lsqdb::LsqDB, saveqr::Dict, Γ; tol_fact = 1.2, verbose=true,
                E0 = nothing, Vref = OneBody(E0))
   R = saveqr["R"]
   Y = saveqr["Y"]
   Nb = size(R, 1)
   y = (Y' * saveqr["Q"])[1:Nb]
   _rmse = norm(Y - saveqr["Q"] * y)
   creg = IPFitting.Lsq.reglsq(; Γ = Γ, R=R, y=y,
                            η0 = _rmse, τ = tol_fact * _rmse,
                            verbose=verbose)
   verbose && @info("|creg| = $(norm(creg)); |Γ creg| = $(norm(Γ*creg))")
   return JuLIP.MLIPs.SumIP(Vref,
                            JuLIP.MLIPs.combine(lsqdb.basis, creg))
end



end

#############
# κ, p_1, int_order = 0.0, 0.0, 0.0
#    if (solver[1] == :qr) || (solver == :qr)
#       verbose && @info("solve $(size(Ψ)) LSQ system using QR factorisation")
#       qrΨ = qr!(Ψ)
#       κ = cond(qrΨ.R)
#       GC.gc();
#       verbose && @info("cond(R) = $(cond(qrΨ.R))")
#       c = qrΨ \ Y
#       rel_rms = norm(qrΨ.Q * (qrΨ.R * c) - Y) / norm(Y)

#       if saveqr isa Dict
#          saveqr["Q"] = qrΨ.Q
#          saveqr["R"] = qrΨ.R
#          saveqr["Y"] = Y
#       end

#       qrΨ = nothing

#    elseif solver[1] == :svd
#       verbose && @info("solve $(size(Ψ)) LSQ system using SVD factorisation")
#       ndiscard = solver[2]
#       F = svd(Ψ)
#       c = F.V[:,1:(end-ndiscard)] * (Diagonal(F.S[1:(end-ndiscard)]) \ (F.U' * Y)[1:(end-ndiscard)])
#       rel_rms = norm(Ψ * c - Y) / norm(Y)

#    elseif solver[1] == :rrqr
#       verbose && @info("solve $(size(Ψ)) LSQ system using Rank-Revealing QR factorisation")
#       qrΨ = pqrfact(Ψ, rtol=solver[2])
#       verbose && @info("cond(R) = $(cond(qrΨ.R))")
#       c = qrΨ \ Y
#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :rid
#       r = solver[2]
#       verbose && @info("solve $(size(Ψ)) LSQ system using Ridge Regression [r = $(r)] ")
#       qrΨ = qr!(Ψ)
#       Nb = size(qrΨ.R, 1)
#       y = (Y' * Matrix(qrΨ.Q))[1:Nb]
#       η0 = norm(Y - Matrix(qrΨ.Q) * y)

#       τ = r * η0
#       Γ = Matrix(I, Nb, Nb)

#       c = reglsq(Γ = Γ, R = Matrix(qrΨ.R), y=y, τ= τ, η0 = η0 );

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :brr
#       BRR = pyimport("sklearn.linear_model")["BayesianRidge"]

#       @info("Using BRR Regression")

#       clf = BRR(normalize=true, compute_score=true)
#       clf.fit(Ψ, Y)

#       c = clf.coef_
#       global alpha = clf.alpha_
#       global beta = clf.lambda_
#       global score = clf.scores_[end]

#       @info("Score: $(score)")

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :ard
#       tol = solver[2]
#       threshold_lambda = solver[3]

#       ARD = pyimport("sklearn.linear_model")["ARDRegression"]
#       @info("Using ARD Regression")
#       @info("Tolerance: $(tol), Threshold lambda: $(threshold_lambda)")

#       clf = ARD(threshold_lambda = threshold_lambda, tol=tol, normalize=true, compute_score=true)
#       clf.fit(Ψ, Y)

#       c = clf.coef_
#       global alpha = clf.alpha_
#       global beta = clf.lambda_
#       global score = clf.scores_[end]

#       zero_ind = findall(x -> x == 0.0, c)
#       non_zero_ind = findall(x -> x != 0.0, c)

#       @info("Fit complete: keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(c), digits=2)*100)%)")
#       @info("Score: $(score)")

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :brr_lap
#       BRR = pyimport("sklearn.linear_model")["BayesianRidge"]

#       rlap_scal = solver[2]
#       #a2b = solver[2][2]
#       @info("rlap_scal=$(rlap_scal)")#, a2b=$(a2b)")

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)#; a2b = a2b)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       clf = BRR(normalize=true)
#       clf.fit(Ψ, Y)

#       creg = clf.coef_
#       score = clf.scores_[end]

#       c = D_inv * creg

#       rel_rms = norm(Ψ * creg - Y) / norm(Y)
#    elseif solver[1] == :lap
#       rscal = solver[2][1]
#       r = solver[2][2]
#       verbose && @info("solve $(size(Ψ)) LSQ system using Laplacian Regularisation [r = $(r)] ")

#       qrΨ = qr!(Ψ)
#       Nb = size(qrΨ.R, 1)
#       y = (Y' * Matrix(qrΨ.Q))[1:Nb]
#       η0 = norm(Y - Matrix(qrΨ.Q) * y)

#       τ = r * η0

#       s = ACE.scaling(db.basis.BB[2], rscal)
#       l = append!(ones(length(db.basis.BB[1])), s)

#       Γ = Diagonal(l)

#       c = reglsq(Γ = Γ, R = Matrix(qrΨ.R), y=y, τ= τ, η0 = η0 );

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :rrqr_lap
#       r_tol = solver[2]
#       rlap = solver[3]

#       s = ACE.scaling(db.basis.BB[2], 2)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(rlap .* l)

#       Ψreg = vcat(Ψ, Γ)

#       qrΨreg = pqrfact(Ψreg, rtol=r_tol)
#       Yreg = vcat(Y, zeros(length(l)))

#       c = qrΨreg \ Yreg

#       rel_rms = norm(qrΨreg * c - Yreg) / norm(Yreg)
#    elseif solver[1] == :rrqr_lap_rescale
#       r_tol = solver[2]
#       rlap = solver[3]

#       s = ACE.scaling(db.basis.BB[2], 2)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(rlap .* l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       qrΨreg = pqrfact(Ψreg, rtol=r_tol)
#       c = D_inv * (qrΨreg \  Y)
#       rel_rms = norm(Ψreg * c - Y) / norm(Y)
#    elseif solver[1] == :lap_elastic_net_rel
#       rlap_scal = solver[2][1]
#       rtol = solver[2][2]
#       etol = solver[2][3]

#       @info("rlap_scal = $(rlap_scal), rrqr_tol=$(rtol), e_tol = $(etol)")

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       non_zero_ind = [j for (j,i) in enumerate(Ψ[1,:]) if sum(i) != 0]
#       zero_ind = [j for (j,i) in enumerate(Ψ[1,:]) if sum(i) == 0]

#       Ψreg_red = Ψreg[:, setdiff(1:end, zero_ind)]

#       qrΨ = pqrfact(Ψreg_red, rtol=rtol)
#       cred = qrΨ \ Y

#       cred_big = zeros(length(Ψreg[1,:]))

#       for (i,k) in enumerate(non_zero_ind)
#         cred_big[k] = cred[i]
#       end

#       c = D_inv * cred_big

#       rel_rms0 = norm(Ψ * c - Y) / norm(Y)

#       @info("rel_rms0 = $(rel_rms0)")
#       @info("etol*rel_rms0 = $(etol*rel_rms0)")

#       function _f(Ψreg, Y, α; etol=1e-5, rtol=1e-9, return_solution=false)
#           cv = glmnet(Ψreg, Y, alpha=α)
#           theta = cv.betas[:, end]

#           non_zero_ind = findall(x -> x != 0.0, theta)
#           zero_ind = findall(x -> x == 0.0, theta)

#           Ψreg_red = Ψreg[:, setdiff(1:end, zero_ind)]

#           qrΨ = pqrfact(Ψreg_red, rtol=rtol)
#           cred = qrΨ \ Y

#           cred_big = zeros(length(Ψreg[1,:]))

#           for (i,k) in enumerate(non_zero_ind)
#               cred_big[k] = cred[i]
#           end

#           c = D_inv * cred_big

#           rel_rms = norm(Ψ * c - Y) / norm(Y)
#           @info("rel_rms=$(rel_rms) and α=$(α), keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(theta), digits=2)*100)%)")

#           if return_solution
#               return c
#           else
#               return rel_rms - etol
#           end
#       end

#       α = find_zero(α -> _f(Ψreg, Y, α, etol=etol*rel_rms0, rtol=rtol), (0,1), Roots.Bisection(), xatol=1E-6)#atol=etol/2) #atol=0.5
#       @info("α found! α=$(α)")

#       c = _f(Ψreg, Y, α, return_solution=true)
#       ##
#       s_1 = ACE.scaling(db.basis.BB[2], 1)
#       p_1 = append!(ones(length(db.basis.BB[1])), s_1)
#       #int_order = db.basis.BB[2].pibasis.inner[1].orders
#       ##
#    elseif solver[1] == :lap_rrqr
#       rlap_scal = solver[2][1]
#       rtol = solver[2][2]
#       reduce = false
#       if length(solver[2]) == 3 && solver[2][3] == true
#          reduce = true
#       end

#       @info("rlap_scal = $(rlap_scal), rrqr_tol=$(rtol)")

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       if reduce
#          Ψreg = convert.(Float32, Ψreg)
#          Y = convert.(Float32, Y)
#       end

#       qrΨ = pqrfact!(Ψreg, rtol=rtol)
#       cred = qrΨ \ Y

#       c = D_inv * cred

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#       ##
#       s_1 = ACE.scaling(db.basis.BB[2], 1)
#       p_1 = append!(ones(length(db.basis.BB[1])), s_1)
#       #int_order = db.basis.BB[2].pibasis.inner[1].orders
#       ##
#    elseif solver[1] == :lap_tik
#       rlap_scal = solver[2][1]
#       r_tik = solver[2][2]

#       @info("rlap_scal = $(rlap_scal), r_tik=$(r_tik)")

#       non_zero_ind = [j for (j,i) in enumerate(Ψ[1,:]) if sum(i) != 0]
#       zero_ind = [j for (j,i) in enumerate(Ψ[1,:]) if sum(i) == 0]

#       nbasis_full = length(Ψ[1,:])
#       Ψ = Ψ[:, setdiff(1:end, zero_ind)]
#       nbasis = length(Ψ[1,:])
#       nobs = length(Ψ[:,1])

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l[:, setdiff(1:end, zero_ind)])

#       D_inv = pinv(Γ)
#       Ψred = similar(Ψ)
#       mul!(Ψred, Ψ, D_inv)
#       Ψ = nothing

#       Ψreg = zeros(nobs + nbasis, nbasis)

#       for (i,row) in enumerate(eachrow(Ψred))
#          Ψreg[i,:] = row
#       end

#       for i in 1:nbasis
#          Ψreg[nobs+i,i] = r_tik
#       end
#       Ψred = nothing

#       Yreg = vcat(Y, zeros(nbasis))

#       GC.gc()
#       qrΨreg = qr!(Ψreg)
#       cred = qrΨreg \ Yreg

#       cred_small = D_inv * cred
#       c = zeros(nbasis_full)

#       for (i,k) in enumerate(non_zero_ind)
#          c[k] = cred_small[i]
#       end

#       rel_rms = 0.0
#       ##
#       s_1 = ACE.scaling(db.basis.BB[2], 1)
#       p_1 = append!(ones(length(db.basis.BB[1])), s_1)
#       #int_order = db.basis.BB[2].pibasis.inner[1].orders
#       ##
#    elseif solver[1] == :elastic_net_lsqr
#       damp = solver[2][1]
#       atol = solver[2][2]
#       α = solver[2][3]

#       @info("damp=$(damp), lsqr_atol=$(atol), alpha=$(α) | ELNET LSQR")

#       cv = glmnet(Ψ, Y, alpha=α)
#       theta = cv.betas[:,end]

#       non_zero_ind = findall(x -> x != 0.0, theta)
#       zero_ind = findall(x -> x == 0.0, theta)

#       @info("keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(theta), digits=2)*100)%)")

#       Ψred = Ψ[:, setdiff(1:end, zero_ind)]

#       cred, lsqrinfo = lsqr(Ψred, Y, damp=damp, atol=atol, log=true)
#       println(lsqrinfo)

#       c = zeros(length(Ψ[1,:]))

#       for (i,k) in enumerate(non_zero_ind)
#         c[k] = cred[i]
#       end

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#       #int_order = db.basis.BB[2].pibasis.inner[1].orders
#       ##
#    elseif solver[1] == :elastic_net_lap
#       α = solver[2][1]
#       rlap_scal = solver[2][2]
#       r = solver[2][3]
#       cv = glmnet(Ψ, Y, alpha=α)
#       theta = cv.betas[:,end]

#       non_zero_ind = findall(x -> x != 0.0, theta)
#       zero_ind = findall(x -> x == 0.0, theta)

#       @info("keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(theta), digits=2)*100)%)")

#       Ψred = Ψ[:, setdiff(1:end, zero_ind)]

#       qrΨred = qr!(Ψred)
#       Nb = size(qrΨred.R, 1)
#       y = (Y' * Matrix(qrΨred.Q))[1:Nb]
#       η0 = norm(Y - Matrix(qrΨred.Q) * y)

#       τ = r * η0

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       lred = [l[i] for i in non_zero_ind]
#       Γ = Diagonal(lred)

#       cred = reglsq(Γ = Γ, R = Matrix(qrΨred.R), y=y, τ= τ, η0 = η0 );

#       c = zeros(length(Ψ[1,:]))

#       for (i,k) in enumerate(non_zero_ind)
#         c[k] = cred[i]
#       end

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :lap_elastic_net_lap
#       α = solver[2][1]
#       rlap_scal = solver[2][2]
#       r = solver[2][3]

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       cv = glmnet(Ψreg, Y, alpha=α)
#       theta = cv.betas[:,end]

#       non_zero_ind = findall(x -> x != 0.0, theta)
#       zero_ind = findall(x -> x == 0.0, theta)

#       Ψred = Ψreg[:, setdiff(1:end, zero_ind)]

#       qrΨred = qr!(Ψred)
#       Nb = size(qrΨred.R, 1)
#       y = (Y' * Matrix(qrΨred.Q))[1:Nb]
#       η0 = norm(Y - Matrix(qrΨred.Q) * y)

#       τ = r * η0

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       lred = [l[i] for i in non_zero_ind]
#       Γ = Diagonal(lred)

#       cred = reglsq(Γ = Γ, R = Matrix(qrΨred.R), y=y, τ= τ, η0 = η0 );

#       cred_big = zeros(length(Ψ[1,:]))

#       for (i,k) in enumerate(non_zero_ind)
#         cred_big[k] = cred[i]
#       end

#       c = D_inv * cred_big

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :gd_reg
#       λ = solver[2][1]
#       rlap_scal = solver[2][2]

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       function error(c, y, Ψ, λ)
#          ŷ = Ψ*c
#          loss = mean((y .- ŷ).^2) + λ*norm(c)
#          return loss
#      end

#       res = optimize(c -> error(c, Y, Ψreg, λ), zeros(length(Ψ[1,:])), GradientDescent())

#       creg = Optim.minimizer(res)

#       c = D_inv * creg

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :lbfgs_reg
#       λ = solver[2][1]
#       rlap_scal = solver[2][2]

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       function _error(c, y, Ψ, λ)
#          ŷ = Ψ*c
#          loss = mean((y .- ŷ).^2) + λ*norm(c)
#          return loss
#       end

#       res = optimize(c -> _error(c, Y, Ψreg, λ), zeros(length(Ψ[1,:])), LBFGS())

#       creg = Optim.minimizer(res)

#       c = D_inv * creg
#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :itlsq
#       damp = solver[2][1]
#       rlap_scal = solver[2][2]
#       atol = solver[2][3]
#       a2b = solver[2][4]
#       if length(solver[2]) == 5
#          maxiter, c_init = solver[2][5]
#          @info("Using a given approximate solution c, maxiter=$(maxiter)")
#       else
#          c_init = zeros(length(db.Ψ[1,:]))
#          maxiter=100000
#       end
#       @info("damp=$(damp), rlap_scal=$(rlap_scal), lsqr_atol=$(atol), a2b=$(a2b)")

#       s = ACE.scaling(db.basis.BB[2], rlap_scal; a2b = a2b)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       creg, lsqrinfo = lsqr!(c_init, Ψ, Y, damp=damp, atol=atol, maxiter=maxiter, log=true)
#       println(lsqrinfo)

#       c = D_inv * creg

#       rel_rms = norm(Ψ * creg - Y) / norm(Y)
#    elseif solver[1] == :itlsq_lap2b
#       damp = solver[2][1]
#       rlap_scal = solver[2][2]
#       atol = solver[2][3]
#       @info("damp=$(damp), rlap_scal=$(rlap_scal), lsqr_atol=$(atol) | scaling 2B ONLY")

#       s = ACE.scaling(db.basis.BB[1], rlap_scal)#; a2b = a2b)
#       l = append!(s, ones(length(db.basis.BB[2])))
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       creg, lsqrinfo = lsqr(Ψ, Y, damp=damp, atol=atol, log=true)
#       println(lsqrinfo)

#       c = D_inv * creg

#       rel_rms = norm(Ψ * creg - Y) / norm(Y)
#    elseif solver[1] == :lap2b_rid
#       r = solver[2][1]
#       rlap_scal = solver[2][2]
#       @info("r=$(r), rlap_scal=$(rlap_scal) | scaling 2B ONLY")

#       s = ACE.scaling(db.basis.BB[1], rlap_scal)#; a2b = a2b)
#       l = append!(s, ones(length(db.basis.BB[2])))
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       qrΨ = qr!(Ψ)
#       Nb = size(qrΨ.R, 1)
#       y = (Y' * Matrix(qrΨ.Q))[1:Nb]
#       η0 = norm(Y - Matrix(qrΨ.Q) * y)

#       τ = r * η0
#       Γ = Matrix(I, length(l), length(l))

#       creg = reglsq(Γ = Γ, R = Matrix(qrΨ.R), y=y, τ= τ, η0 = η0 );

#       c = D_inv * creg

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :lsqr
#       damp = solver[2][1]
#       atol = solver[2][2]
#       @info("damp=$(damp), lsqr_atol=$(atol) | SIMPLE LSQR")

#       c, lsqrinfo = lsqr(Ψ, Y, damp=damp, atol=atol, log=true)
#       println(lsqrinfo)

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :itlsq_Nw
#       damp = solver[2][1]
#       rlap_scal = solver[2][2]
#       atol = solver[2][3]
#       Nweight = solver[2][4]
#       @info("damp=$(damp), rlap_scal=$(rlap_scal), lsqr_atol=$(atol)")

#       V2scal = Nweight["2B"] .* ACE.scaling(db.basis.BB[1], rlap_scal)
#       NBscal = ACE.scaling(db.basis.BB[2], rlap_scal)
#       RegDiag = vcat(V2scal, NBscal)
#       Γ = Diagonal(RegDiag)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       creg, lsqrinfo = lsqr(Ψ, Y, damp=damp, atol=atol, log=true)
#       println(lsqrinfo)

#       c = D_inv * creg

#       rel_rms = norm(Ψ * creg - Y) / norm(Y)
#    elseif solver[1] == :elastic_net_itlsq
#       α = solver[2][1]
#       damp = solver[2][2]
#       rlap_scal = solver[2][3]
#       atol = solver[2][4]
#       a2b = solver[2][5]
#       @info("α=$(α), damp=$(damp), rlap_scal=$(rlap_scal), lsqr_atol=$(atol), a2b=$(a2b)")

#       s = ACE.scaling(db.basis.BB[2], rlap_scal; a2b = a2b)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       cv = glmnet(Ψ, Y, alpha=α)
#       theta = cv.betas[:,end]

#       non_zero_ind = findall(x -> x != 0.0, theta)
#       zero_ind = findall(x -> x == 0.0, theta)

#       @info("α=$(α), keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(theta), digits=2)*100)%)")
#       Ψ = Ψ[:, setdiff(1:end, zero_ind)]

#       cred, lsqrinfo = lsqr(Ψ, Y, damp=damp, atol=atol, log=true)
#       println(lsqrinfo)

#       cred_big = zeros(length(l))

#       for (i,k) in enumerate(non_zero_ind)
#         cred_big[k] = cred[i]
#       end

#       c = D_inv * cred_big

#       rel_rms = norm(Ψ * cred - Y) / norm(Y)
#       #rel_rms = 0.0
#       rel_rms = norm(Ψ * cred - Y) / norm(Y)
#    elseif solver[1] == :no2b_itlsq
#       damp = solver[2][1]
#       rlap_scal = solver[2][2]
#       atol = solver[2][3]
#       @info("Not keeping ACE 2B, damp=$(damp), rlap_scal=$(rlap_scal), lsqr_atol=$(atol)")

#       len_pair = length(db.basis.BB[1])
#       len_ace = length(db.basis.BB[2])

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(len_pair), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       mul!(Ψ,Ψ,D_inv)

#       I_ord2 = findall(get_orders(db.basis.BB[2]) .== 1);

#       zero_ind = I_ord2 .+ len_pair
#       non_zero_ind = []
#       for i in 1:(len_pair+len_ace)
#           if i ∉ zero_ind
#               append!(non_zero_ind, i)
#           end
#       end

#       #theta = append!(ones(length(db.basis.BB[1])), cs)

#       #non_zero_ind = findall(x -> x != 0.0, theta)
#       #zero_ind = findall(x -> x == 0.0, theta)

#       @info("keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/(len_pair+len_ace), digits=2)*100)%)")
#       Ψ = Ψ[:, setdiff(1:end, zero_ind)]

#       cred, lsqrinfo = lsqr(Ψ, Y, damp=damp, atol=atol, log=true)
#       println(lsqrinfo)

#       cred_big = zeros(length(l))

#       for (i,k) in enumerate(non_zero_ind)
#         cred_big[k] = cred[i]
#       end

#       c = D_inv * cred_big

#       rel_rms = norm(Ψ * cred - Y) / norm(Y)
#       #rel_rms = 0.0
#       rel_rms = norm(Ψ * cred - Y) / norm(Y)
#    elseif solver[1] == :elastic_net_lap_manual
#       α = solver[2][1]
#       rscal = solver[2][2]
#       lap_pen = solver[2][3]

#       cv = glmnet(Ψ, Y, alpha=α)
#       theta = cv.betas[:,end]

#       non_zero_ind = findall(x -> x != 0.0, theta)
#       zero_ind = findall(x -> x == 0.0, theta)

#       @info("keeping $(length(non_zero_ind)) basis functions ($(round(length(non_zero_ind)/length(theta), digits=2)*100)%)")

#       Ψred = Ψ[:, setdiff(1:end, zero_ind)]

#       s = ACE.scaling(db.basis.BB[2], rscal)
#       l = append!(ones(length(db.basis.BB[1])), s)

#       lred = [l[i] for i in non_zero_ind]
#       Γ = Diagonal(lap_pen .* lred)

#       Ψreg = vcat(Ψred, collect(Γ))
#       Yreg = vcat(Y, zeros(length(lred)))

#       #qrΨreg = qr!(Ψreg)
#       cred = Ψreg \ Yreg

#       c = zeros(length(Ψ[1,:]))

#       for (i,k) in enumerate(non_zero_ind)
#         c[k] = cred[i]
#       end

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    elseif solver[1] == :spgl_lap_rrqr
#       it = solver[2][1]
#       rlap_scal = solver[2][2]
#       rtol = solver[2][3]

#       s = ACE.scaling(db.basis.BB[2], rlap_scal)
#       l = append!(ones(length(db.basis.BB[1])), s)
#       Γ = Diagonal(l)

#       D_inv = pinv(Γ)
#       Ψreg = Ψ * D_inv

#       opts = spgOptions(iterations = it)
#       theta, r, g, info = spgl1(Ψ, Y, sigma = 1e-5, options=opts)

#       non_zero_ind = findall(x -> x != 0.0, theta)
#       zero_ind = findall(x -> x == 0.0, theta)

#       Ψred = Ψreg[:, setdiff(1:end, zero_ind)]

#       qrΨ = pqrfact!(Ψred, rtol=rtol)
#       cred = qrΨ \ Y

#       cred_big = zeros(length(Ψ[1,:]))

#       @show length(cred_big)
#       @show length(zero_ind)

#       for (i,k) in enumerate(non_zero_ind)
#         cred_big[k] = cred[i]
#       end

#       c = D_inv * cred_big

#       rel_rms = norm(Ψ * c - Y) / norm(Y)
#    else
#       error("unknown `solver` in `lsqfit`")
#    end



# @noinline function onb(db::LsqDB;
#                          solver=(:qr, ), verbose=true,
#                          Ibasis = :,
#                          Itrain = :,
#                          E0 = nothing,
#                          Vref = OneBody(E0),
#                          weights = nothing,
#                          regularisers = [],
#                          combineIP = nothing,
#                          kwargs...)
#
#    verbose && @info("assemble lsq system")
#    Ψ, _ = get_lsq_system(db; verbose=verbose, Vref=Vref,
#                              Ibasis=Ibasis, Itrain=Itrain,
#                              weights = weights,
#                              regularisers = regularisers,
#                              kwargs...)
#    @assert solver[1] == :qr
#    verbose && @info("QR-factorize Ψ, size=$(size(Ψ))")
#    qrΨ = qr(Ψ)
#    verbose && @info("cond(R) = $(cond(qrΨ.R))")
#    Rinv = pinv(qrΨ.R)
#    basis = db.basis[Ibasis]
#    onb = []
#    for n = 1:size(Rinv, 2)
#       push!(onb, combineIP(basis, Rinv[:,n]))
#    end
#    return [b for b in onb]
# end
#
# rlap = solver[2][2]
# p = solver[2][3]
#
# scale = true
# if length(solver[2]) == 4 && solver[2][4] == false
#    scale = false
# end
#
# s_p = scaling(db.basis.BB[2], p)
# p_p = append!(ones(length(db.basis.BB[1])), s_p)
# Γ = collect(Diagonal(rlap .* p_p))
#
# if scale == true
#    D_inv = pinv(Γ)
#    Ψreg = Ψ * D_inv
#    Yreg = Y
# else
#    Ψreg = vcat(Ψ, Γ)
#    Yreg = vcat(Y, ones(length(Ψ[1,:])))
# end
#
# if startswith(String(solver[1]), "lap_elastic_net_perc")
#    perc_ob = solver[2][1]
#    function _f(α, perc_ob, Ψreg, Yreg)
#       cv = glmnet(Ψreg, Yreg, alpha=α)
#       theta = cv.betas[:,end]
#
#       zero_ind = findall(x -> x == 0.0, theta)
#       non_zero_ind = findall(x -> x != 0.0, theta)
#
#       perc = round(length(non_zero_ind)/length(theta) * 100.0, digits=2)
#       @info("α=$(α) keeps $(perc)% of total basis functions")
#
#       return length(non_zero_ind)/length(theta) - perc_ob
#    end
#
#    α = find_zero(α -> _f(α, perc_ob, Ψreg, Yreg), (0,1), Roots.Bisection(), atol=0.01)
#    @info("α found! α=$(α)")
# else
#    α = solver[2][1]
# end
#
# cv = glmnet(Ψreg, Yreg, alpha=α)
# theta = cv.betas[:,end]
#
# zero_ind = findall(x -> x == 0.0, theta)
# non_zero_ind = findall(x -> x != 0.0, theta)
#
# perc = round(length(non_zero_ind)/length(db.basis)*100, digits=2)
# @info("Reduced LSQ Problem using $(perc)% of total basis functions [$(length(non_zero_ind))]")
#
# if scale == true
#    Ψred = Ψreg[:, setdiff(1:end, zero_ind)]
# else
#    Ψred = Ψ[:, setdiff(1:end, zero_ind)]
# end
#
# if endswith(String(solver[1]), "rrqr")
#    rtol = solver[3]
#    @info("Performing RRQR [$(rtol)]")
#    qrΨred = pqrfact(Ψred, rtol=rtol)
#    cred = qrΨred \ Y
# else
#    @info("Performing QR decomposition")
#    cred = Ψred \ Y
# end
#
# c_scal = zeros(length(db.basis))
#
# for (i,k) in enumerate(non_zero_ind)
#    c_scal[k] = cred[i]
# end
#
# if scale == true
#    c = D_inv * c_scal
# else
#    c = c_scal
# end
#
# if any(iszero, c[1:length(db.basis.BB[1])])
#    @info("Warning: culled a 2B function!! [TODO]")
# end
#
# rel_rms = norm(Ψ * c - Y) / norm(Y)
#
# ##
# s_1 = scaling(db.basis.BB[2], 1)
# p_1 = append!(ones(length(db.basis.BB[1])), s_1)
# ##

