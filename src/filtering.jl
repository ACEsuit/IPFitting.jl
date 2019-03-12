
module Filtering

using NBodyIPFitting: LsqDB
export filter_basis, filter_configs,
       anyf, allf


function filter_basis(db::LsqDB, args...)
   @assert length(args) > 0
   choose_b = fill(false, length(db.basis))
   for (i, b) in enumerate(db.basis)
      choose_b[i] = all( a(b) for a in args )
   end
   return findall(choose_b)
end


function filter_configs(db::LsqDB, args...)
   @assert length(args) > 0
   choose_cfg = fill(false, length(db.configs))
   for (i, cfg) in enumerate(db.configs)
      choose_cfg[i] = all( a(cfg) for a in args )
   end
   return findall(choose_cfg)
end

anyf(args...) = x -> any( f(x) for f in args )
allf(args...) = x -> all( f(x) for f in args )


# TODO: create a random filter 

end
