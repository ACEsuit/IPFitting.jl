
module Tools

using ProgressMeter, Base.Threads

"""
`tfor(f, rg; verbose=true, msg="tfor", costs = ones(Int, length(rg)))`

Multi-threaded for loop. At each iteration the function f(n) is executed,
where n loops over `rg`.
"""
function tfor(f, rg; verbose=true, msg="tfor", costs = ones(Int, length(rg)),
                     maxnthreads=nthreads())
   p = Progress(sum(costs))
   p_ctr = 0
   t0 = time_ns()
   nthr = max(1, min(nthreads(), maxnthreads))
   if nthr == 1
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
         ProgressMeter.update!(p, 0)  # not sure this is useful/needed
         p_lock = SpinLock()
      end
      # sort the tasks by cost: do the expensive ones first
      Isort = sortperm(costs, rev=true)
      # remember what the last job was
      last_job = 0
      last_job_lock = SpinLock()
      # start a simple loop over nthreads() just to split this
      # into parallel "tasks"
      rgidx = Vector{Int}(undef, nthreads())
      @threads for i = 1:nthreads()
         while true
            # acquire a new job index
            tid = threadid()
            lock(last_job_lock)
            last_job += 1
            if last_job > length(Isort)
               unlock(last_job_lock)
               break # break the while loop and then wait for the
                     # other threads to finish
            end
            rgidx[tid] = Isort[last_job]
            unlock(last_job_lock)
            # do the job
            f(rg[rgidx[tid]])
            # submit progress
            if verbose
               lock(p_lock)
               p_ctr += costs[rgidx[tid]]
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



function pfor(db, f, rg; verbose=true, msg="tfor", costs = ones(Int, length(rg)),
                     maxprocs=nprocs()-1)
   p = Progress(sum(costs))
   p_ctr = 0
   t0 = time_ns()
   nprcs = max(1, min(nprocs() - 1, maxprocs))
   #TODO: serial will not work for now
   if nprcs == 1
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
         @info("$msg with $(nprcs) processors")
         ProgressMeter.update!(p, 0)  # not sure this is useful/needed
         p_lock = SpinLock()
      end
      # sort the tasks by cost: do the expensive ones first
      Isort = sortperm(costs, rev=true)

      # remember what the last job was
      last_job = 1
      not_done = true
      #stores the current work the processors are doing
      futures = Vector{Future}(undef, nprcs)

      # loop until all jobs are done
      while not_done
         # spawn the first job for all workers
         for i in 1:nprcs
            tmpf = @spawnat (i+1) f(rg[Isort[last_job]])
            futures[i] = tmpf
            # submit progress
            if verbose
               p_ctr += costs[Isort[last_job]]
               ProgressMeter.update!(p, p_ctr)
            end
            last_job += 1
            # check if we are done
            if last_job > length(Isort)
               not_done = false
               break
            end
         end

         # loop through the workers, if they are done, record that value in the matrix and give them
         # another job
         for i in 1:nprcs
            if(isready(futures[i]))
               (irows,vec_lsqrow) = fetch(futures[i])
               db.Ψ[irows, :] = vec_lsqrow
               tmpf = @spawnat (i+1) f(rg[Isort[last_job]])
               futures[i] = tmpf
               # submit progress
               if verbose
                  p_ctr += costs[Isort[last_job]]
                  ProgressMeter.update!(p, p_ctr)
               end
               last_job += 1
               if last_job > length(Isort)
                  not_done = false
                  break 
               end
            end
         end
      end

      # Once done we go over the features to get the last ones
      # so far this is naive and will double update some of them TODO: improve this
      for i in 1:nprcs
         (irows,vec_lsqrow) = fetch(futures[i])
         db.Ψ[irows, :] = vec_lsqrow
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
