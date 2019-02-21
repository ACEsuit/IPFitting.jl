module Weights

import JuLIP, NBodyIPFitting
using NBodyIPFitting: LsqDB, Dat, configname
using JuLIP: forces, Atoms, positions

function hess_weights!(cn::String, db::LsqDB; h = :auto, wE = 0.0, rscal = 3, verbose=false)

   # temporarily assume h is a floating point number
   # => this needs to be automated
   @assert (h isa AbstractFloat)

   # loop through the data_groups
   for (key, dg) in db.data_groups
      # if the data_group has the correct configname
      if cn == configname(key)
         # then loop through the dats
         for d in dg
            if forces(d) == nothing
               verbose && @warn("""hess_weights!: a training configuration does
                                  not contain energy and forces => ignore""")
               d.info["W"] = Dict( "F" => zeros(3*length(d)),
                                   "E" => 0.0,
                                   "V" => zeros(6) )
               continue
            end
            # get configuration info, R vectors
            at = Atoms(d)
            X = positions(at)
            R = [ JuLIP.project_min(at, x - X[1])  for x in X ]

            # if h == :auto => COMPUTE h <<<<<<<<

            # now fix the scaling for the force weights
            r = norm.(R)
            r3 = (ones(3) * r')[:]
            wF = r3.^rscal / h
            wF[1:3] = 0.0   # this sets the weight for atom 1 to zero
                            # it contains no information...
            # assign this weight back to `d::Dat`
            d.info["W"] = Dict( "E" => wE,
                                "F" => wF,
                                "V" => zeros(6) )
         end
      end
   end

   return nothing
end


end # module
