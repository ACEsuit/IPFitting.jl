using Combinatorics, StaticArrays

include("misc.jl")

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

   Output = [];

   open(filename, "a") do f
      for j=1:deg
         write(f, pref)
         write(f, "_$j = ")
         write(f, "@SVector [")
         push!(Output, pref*"_$j")
         for i=1:(L-1)
            k = Ind[i,j];
            write(f, "x[$k], ")
         end
         k = Ind[L,j];
         write(f, "x[$k] ] \n" )
      end
   end
   return Output
end

# function monomial_2_file(Perms,Permsref,filename,pref)
filename = "test2.jl";
pref = "PV1";
x = @SVector [3,1,0,0,0,0]
xref = @SVector [1,1,0,0,0,0]
Permsref = xref
Perms = x

   # nb of non-zero coefficients
   perm_deg = sum(Permsref)
   nb_diff_coef = length(unique(Perms))-1
   Ind = perm_2_indice(unique(simplex_permutations(Permsref)))
   Names = indice_2_file(Ind,filename,pref)
   Decreasing_coef = sort(Perms,rev=true)
   if (perm_deg==1)
      open(filename, "a") do f
         exponent = Decreasing_coef[1];
         write(f, pref)
         write(f, " = ")
         write(f, "sum(x$exponent)\n")
      end
   elseif (maximum(Perms)==1)&(perm_deg==2)
      open(filename, "a") do f
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*","Names[2]")\n")
      end
   elseif (maximum(Perms)==1)&(perm_deg==3)
      open(filename, "a") do f
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*".*"Names[2]","*Names[3]")\n")
      end
   elseif (maximum(Perms)==1)&(perm_deg==4)
      open(filename, "a") do f
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*".*"Names[2]","*Names[3]".*"*Names[4]")\n")
      end
   elseif (maximum(Perms)==1)&(perm_deg==5)
      open(filename, "a") do f
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*".*"Names[2]","*Names[3]".*"*Names[4]".*"*Names[5]")\n")
      end
   elseif (Decreasing_coef[1]==2)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_2,"*Names[2]*")+dot("*Names[1]","Names[2]"_2)\n")
      end
   elseif (Decreasing_coef[1]==3)&(Decreasing_coef[2]==1)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_3,"*Names[2]*")+dot("*Names[1]","Names[2]"_3)\n")
      end
   else
      error("not implemented yet")
   end
# end

x = @SVector [2,1,0,0,0,0]
xref = @SVector [1,1,0,0,0,0]

# y = unique(simplex_permutations(x))
# length(y)
# sum(y[1])
# M = length(y[1])
# DD = find(y[1])
# Ind = perm_2_indice(y)
#
#
# indice_2_file(Ind,"test.jl","I4")

monomial_2_file(x,xref,"test2.jl","PV1")
