
import Base.iterate

struct ObservationsIterator
   configs::Vector{Dat}
end

observations(db::LsqDB) = observations(db.configs)
observations(configs::Vector{Dat}) = ObservationsIterator(configs)

iterate(iter::ObservationsIterator) = iterate(iter, 1)

function iterate(iter::ObservationsIterator, i::Integer)
   if i > length(iter.configs)
      return nothing
   end
   # the observation keys for the current config
   okeys = collect(keys(iter.configs[i].D))
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
   d = iter.configs[state.i]
   return (obskey, d, state.i),  # returned to user
          (i=state.i, okeys=state.okeys, ikey=state.ikey+1)   # next state
end


tfor_observations(db::LsqDB, callback; kwargs...) =
      observations(db.configs, callback; kwargs...)

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
