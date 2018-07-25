
module Tools

function tfor(f, rg; verbose=true, msg="")
   if nthreads() == 1
      verbose && println("$msg in serial")
      dt = verbose ? 1.0 : Inf
      tic()
      @showprogress dt for n in rg
         f(n)
      end
      verbose && toc()
   else
      if verbose
         println("$msg with $(nthreads()) threads")
         p = Progress(length(rg))
         p_ctr = 0
         p_lock = SpinLock()
      end
      tic()
      @threads for n in rg
         f(n)
         if verbose
            lock(p_lock)
            p_ctr += 1
            ProgressMeter.update!(p, p_ctr)
            unlock(p_lock)
         end
      end
      verbose && toc()
   end
   return nothing
end

end
