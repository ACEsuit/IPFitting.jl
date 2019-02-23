
module Tools

using ProgressMeter, Base.Threads

function tfor(f, rg; verbose=true, msg="tfor", costs = ones(Int, length(rg)))
   p = Progress(sum(costs))
   p_ctr = 0
   if nthreads() == 1
      verbose && println("$msg in serial")
      dt = verbose ? 1.0 : Inf
      t0 = time_ns()
      for (n, c) in zip(rg, costs)
         f(n)
         if verbose
            p_ctr += c
            ProgressMeter.update!(p, p_ctr)
         end
      end
      t1 = time_ns()
      verbose && @info("Elapsed: $(round((t1-t0)*1e-9, digits=1))s")
   else
      if verbose
         @info("$msg with $(nthreads()) threads")
         p_lock = SpinLock()
      end
      # tic()
      @threads for i = 1:length(rg)
         f(rg[i])
         if verbose
            lock(p_lock)
            p_ctr += costs[i]
            ProgressMeter.update!(p, p_ctr)
            unlock(p_lock)
         end
      end
      # verbose && toc()
   end
   return nothing
end

decode(D::Dict) = convert(Val(Symbol(D["__id__"])), D)



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
