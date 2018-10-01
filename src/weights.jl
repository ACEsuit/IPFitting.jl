module Weights

using NBodyIPFitting: LsqDB, Dat

function hess_weights!(cn::String, db::LsqDB; h = :auto, wE = 0.0, rscal = 3)

   # loop through the data_groups
   for dg in db.data_groups
      # if the data_group has the correct configname
      if cn == configname(key(dg))
         # then loop through the dats
         for d in dg
            if forces(d) == nothing
               warn("""hess_weights!: a training configuration does not contain
                       energy and forces => ignore""")
               return zeros(length(d))
            end
            # get configuration info, R vectors
            at = Atoms(d)
            X = positions(at)
            R = [ JuLIP.project_min(at, x - X[1])  for x in X ]

            # COMPUTE h <<<<<<<<

            # now fix the scaling for the force weights
            r = norm.(R)
            r3 = (ones(3) * r')[:]
            w = [ wE;
                  r3.^rscal / h ]
            w[2:4] = 0.0    # this sets the weight for atom 1 to zero
                            # it contains no information...
            # assign this weight back to `d::Dat`
            d.D["W"] = w
         end
      end
   end

   return nothing
end


end # module
