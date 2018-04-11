"""
return a list of all the unique permutations of a vector alpha
"""
function uniqueperms(alpha::Vector{Int64})
   size = length(alpha)
   sort!(alpha)
   isFinished = 0
   perms = Vector{Int64}[]
   it = 0
   while isFinished == 0
      push!(perms, copy(alpha))
      k = -1
      for i = (size-1):-1:1
         if alpha[i] < alpha[i+1]
            k = i
            break
         end
      end
      if k == -1
         isFinished = 1
         break
      else
         ceilIndex = k + 1
         for i = k + 2:size
            if alpha[ceilIndex] > alpha[i] && alpha[i] > alpha[k]
               ceilIndex = i
            end
         end
         t = alpha[ceilIndex]
         alpha[ceilIndex] = alpha[k]
         alpha[k] = t
      end
      temp = zeros(Int64,size-k)
      for i = 1:(size-k)
         temp[i] = alpha[i + k]
      end
      sort!(temp)
      for i = 1:(size-k)
         alpha[i + k] = temp[i]
      end
   end
   return perms
end

"""
`psym_monomial(alpha, dict, var) -> Expr, Function`

Function to construct the symbol expression and a wrapped function for the
monomial symmetric polynomial corresponding to the vector alpha. We replace
`x[i]^1` by `x[i]`, which may be pointless/slow. Do not replace `x[i]^0` by `1`,
`:(1)` is not treated as an expression. Check that the letter used for the
variable is not used elsewhere (e.g x in exp) as the variable is automatically
replaced.

### Example
```
ex, f = psym_monomial([0,1,1], ["y^0","y^1","y^2"], "y")
f([1., 2., 3.])
```
"""
function psym_monomial(alpha, dict, sym; simplify = true)
   dict = ["1"; dict]
   dim = length(alpha)
   ex = ""
	for p in uniqueperms(alpha)
      ext = ""
      for i = 1:dim
			ext = "$ext*" * replace(dict[p[i]+1], sym, "$(sym)$i")
		end
      ex = ex * "+$(ext[2:end])"
	end
	ex = parse(ex[2:end])

   if simplify
      exs = Calculus.simplify(ex)
   else
      exs = ex
   end
   s = Symbol(sym)
   f = eval(:($s -> $(ind2vec(exs, dim, sym))))
   df = x -> ForwardDiff.gradient(f, x)
	return ex, f, df
end

psym_monomial(alpha, t::Tuple; kwargs...) =
      psym_monomial(alpha, t[1], t[2]; kwargs...)


"""
`psym_polys(dim, dict, sym):`

collect all permutation invariant polynomials up to degree `deg`
in dimension dim. This assembles up to

### Example
```
exs, fs = PermPolys(3, ["y^0","y^1","y^2"], "y")
fs[1]([1., 2., 3.])
```
"""
function psym_polys(dim::Integer, dict, sym; simplify = true)
	polys_ex = Expr[]
	polys_f = Function[]
	polys_df = Function[]
	for i in 1:length(dict)
		for m = 1:dim, alpha in collect(partitions(i, m))
         append!(alpha, zeros(Int, dim - length(alpha)))
         mex, mf, mdf = psym_monomial(alpha, dict, sym; simplify=simplify)
         push!(polys_ex, mex)
         push!(polys_f, mf)
         push!(polys_df, mdf)
      end
	end
	return polys_ex, polys_f, polys_df
end


"""
which dimensionality corresponds to a body-order
"""
nbody_dim(bo::Integer) = (bo * (bo-1)) ÷ 2

"""
`psym_polys_nbody(bo::Integer, dict, sym)`

* `bo` : body order
* `dict` : 1D dictionary
* `sym` : symbol used in the dictionary
"""
function psym_polys_nbody(N::Integer, dict, sym; simplify = true)
	polys_ex = Expr[]
	polys_f = Function[]
	polys_df = Function[]
   # get the lower and upper dimensionality for genuine N-body terms
   dim_lo = nbody_dim(N-1)+1
   dim_hi = nbody_dim(N)
   if length(dict) < dim_lo
      warn("the length of the dictionary is too short for $N-body terms")
   end
	for i in dim_lo:length(dict)
		for m = dim_lo:dim_hi, alpha in collect(partitions(i, m))
         append!(alpha, zeros(Int, dim_hi - length(alpha)))
         mex, mf, mdf = psym_monomial(alpha, dict, sym; simplify=simplify)
         push!(polys_ex, mex)
         push!(polys_f, mf)
         push!(polys_df, mdf)
      end
	end
	return polys_ex, polys_f, polys_df
end
