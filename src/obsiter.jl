
import Base.iterate
using IPFitting.Tools: tfor
using Base.Threads: SpinLock, nthreads

struct ObservationsIterator
   configs::Vector{Dat}   # remember the database
   Icfg::Vector{Int}      # vector of indices over which we iterate
end

"""
```
observations(db::LsqDB, Icfg=:)
observations(configs::Vector{Dat}, Icfg=:)
```

Creates an iterator over all observations stored in the
list of configurations. I.e., the resulting loop will iterate
over all configurations and within each configuration over all
observations stored in that configuration. E.g.,
```
for (okey, cfg, icfg) in observations(configs)
   # okey -> String key
   # cfg :: Dat
   # icfg :: Int, the index of cfg in configs
end
```
"""
observations(db::LsqDB, args...) = observations(db.configs, args...)
observations(configs::Vector{Dat}) = observations(configs, 1:length(configs))
observations(configs::Vector{Dat}, Icfg) =
     ObservationsIterator(configs, collect(Icfg))
observations(configs::Vector{Dat}, ::Colon) =
     observations(configs, 1:length(configs))

iterate(iter::ObservationsIterator) = iterate(iter, 1)

function iterate(iter::ObservationsIterator, i::Integer)
   if i > length(iter.Icfg)
      return nothing
   end
   # the observation keys for the current config
   okeys = collect(keys(iter.configs[iter.Icfg[i]].D))
   # get the next observation
   return iterate(iter, (i=i, okeys=okeys, ikey=1))
end

function iterate(iter::ObservationsIterator, state::NamedTuple)
   # if we got through all observations, then advance to the next config
   if state.ikey > length(state.okeys)
      return iterate(iter, state.i+1)
   end
   # if there is another observation left...
   obskey = state.okeys[state.ikey]
   icfg = iter.Icfg[state.i]
   d = iter.configs[icfg]
   return (obskey, d, icfg),  # returned to user
          (i=state.i, okeys=state.okeys, ikey=state.ikey+1)   # next state
end


tfor_observations(db::LsqDB, callback; kwargs...) =
      tfor_observations(db.configs, callback; kwargs...)

"""
```
tfor_observations(configs::Vector{Dat}, callback; kwargs...)
tfor_observations(db::LsqDB, callback; kwargs...)
```

Create a multi-threaded for loop over all observations in the db of
list of configs. The callback is expected to be of the form
```
callback(n, obskey::AbstractString, cfg::Dat, lock::SpinLock)
```

### Kwargs
* `verbose=true`
* `msg = "Loop over observations"`
"""
function tfor_observations(configs::Vector{Dat}, callback;
                           verbose=true,
                           msg = "Loop over observations",
                           maxnthreads=nthreads())

   # collect a complete list of observations
   idats = Int[]; sizehint!(idats, 3*length(configs))
   okeys = String[]; sizehint!(okeys, 3*length(configs))
   for (okey, _, idat) in observations(configs)
      push!(idats, idat)
      push!(okeys, okey)
   end
   # a rough cost estimate
   costs = [ length(configs[i]) for i in idats ]

   # start a threaded for loop
   lck = SpinLock()
   tfor( n -> callback(n, okeys[n], configs[idats[n]], lck),
         1:length(idats),
         verbose=verbose, msg = msg,
         costs = costs,
         maxnthreads=maxnthreads )

   return nothing
end
