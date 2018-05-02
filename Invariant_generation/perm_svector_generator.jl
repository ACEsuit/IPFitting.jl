using Combinatorics, StaticArrays, NBodyIPs

NBody = 4;
NBlengths = NBody*(NBody-1)/2;

function perm_2_indice(Perms)
   L = length(Perms);
   deg = sum(Perms[1]);
   M = length(Perms[1]);

   Indices = Array{Int64, 2}(L,deg);
   # zeros(L,deg);

      for j=1:L
         Indices[j,:] = find(Perms[j])
      end
   return Indices
end




function indice_2_file(Ind,filename,pref)
   L = size(Ind,1);
   deg = size(Ind,2);

   open(filename, "w") do f
      for j=1:deg
         write(f, pref)
         write(f, "_$j = ")
         write(f, "@SVector [")
         for i=1:(L-1)
            k = Ind[i,j];
            write(f, "x[$k], ")
         end
         k = Ind[L,j];
         write(f, "x[$k] ] \n" )
      end
   end
end

function monomial_2_file(Perms,filename,pref)
   if maximum(Perms)==1
      Ind = perm_2_indice(unique(simplex_permutations(Perms)))
      indice_2_file(Ind,filename,pref)
   else

   end
end

x = @SVector [1,1,1,0,0,0]

y = unique(simplex_permutations(x))
length(y)
sum(y[1])
M = length(y[1])
DD = find(y[1])
Ind = perm_2_indice(y)

indice_2_file(Ind,"test.jl","I4")

monomial_2_file(x,"test2.jl","I4")
