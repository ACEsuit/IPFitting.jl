
module Tools

using ProgressMeter, Base.Threads

"""
`tfor(f, rg; verbose=true, msg="tfor", costs = ones(Int, length(rg)))`

Multi-threaded for loop. At each iteration the function f(n) is executed,
where n loops over `rg`.
"""
function tfor(f, rg; verbose=true, msg="tfor", costs = ones(Int, length(rg)))
   p = Progress(sum(costs))
   p_ctr = 0
   t0 = time_ns()
   if nthreads() == 1
      verbose && println("$msg in serial")
      dt = verbose ? 1.0 : Inf
      for (n, c) in zip(rg, costs)
         f(n)
         if verbose
            p_ctr += c
            ProgressMeter.update!(p, p_ctr)
         end
      end
   else
      if verbose
         @info("$msg with $(nthreads()) threads")
         p_lock = SpinLock()
      end
      # sort the tasks by cost: do the expensive ones first
      Isort = sortperm(costs, rev=true)
      # remember what the last job was
      last_job = 0
      last_job_lock = SpinLock()
      # start a simple loop over nthreads() just to split this
      # into parallel "tasks"
      @threads for i = 1:nthreads()
         while true
            # acquire a new job index
            lock(last_job_lock)
            last_job += 1
            if last_job > length(Isort)
               unlock(last_job_lock)
               break # break the while loop and then wait for the
                     # other threads to finish
            end
            local rgidx = Isort[last_job]
            unlock(last_job_lock)
            # do the job
            f(rg[rgidx])
            # submit progress
            if verbose
               lock(p_lock)
               p_ctr += costs[rgidx]
               ProgressMeter.update!(p, p_ctr)
               unlock(p_lock)
            end
         end
      end
   end
   t1 = time_ns()
   verbose && @info("Elapsed: $(round((t1-t0)*1e-9, digits=1))s")
   return nothing
end


function analyse_include_exclude(set, include, exclude)
   if include != nothing && exclude != nothing
      error("only one of `include`, `exclude` may be different from `nothing`")
   end
   if include != nothing
      if !issubset(include, set)
         error("`include` can only contain elements of `set`")
      end
      # do nothing - just keep `include` as is to return
   elseif exclude != nothing
      if !issubset(exclude, set)
         error("`exclude` can only contain config types that are in `set`")
      end
      include = setdiff(set, exclude)
   else
      # both are nothing => keep all config_types
      include = set
   end
   return include
end


macro def(name, definition)
    return quote
        macro $(esc(name))()
            esc($(Expr(:quote, definition)))
        end
    end
end


end
