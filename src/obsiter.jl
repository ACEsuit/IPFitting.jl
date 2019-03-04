
import Base.iterate

struct ObservationsIterator
   configs::Vector{Dat}   # remember the database
   Icfg::Vector{Int}      # vector of indices over which we iterate
end

observations(db::LsqDB, args...) = observations(db.configs, args...)
observations(configs::Vector{Dat}) = observations(configs, 1:length(configs))
observations(configs::Vector{Dat}, Icfg) =
     ObservationsIterator(configs, collect(Icfg))

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

function tfor_observations(configs::Vector{Dat}, callback;
                  verbose=true,
                  msg = "Loop over observations")

   # collect a complete list of observations
   idats = Int[]; sizehint!(idats, 3*length(configs))
   okeys = String[]; sizehint!(okeys, 3*length(configs))
   for (okey, _, idat) in observations(configs)
      push!(idats, idat)
      push!(okeys, okey)
   end

   # start a threaded for loop
   db_lock = SpinLock()
   tfor( n -> callback(okeys[n], configs[idats[n]], db, db_lock),
         1:length(idats),
         verbose=verbose, msg = msg,
         costs = [ length(configs[i]) for i in idats ] )

   return nothing
end
