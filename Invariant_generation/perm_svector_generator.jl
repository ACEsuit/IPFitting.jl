using Combinatorics, StaticArrays, NBodyIPs

# include("misc.jl")



# TODO: code this function in a nicer way
function generate_monomials_irr_sec(NBody,Deg)
    NBlengths = Int(NBody*(NBody-1)/2);
    NB_sec_inv = countlines("Invariant_generation/NBody_$NBody""_deg_$Deg""_invariants.jl");

    Monomials = Array{Int64, 2}(NB_sec_inv,NBlengths);
    Monomials[:] = 0;

    Monomials_simple = Array{Int64, 2}(NB_sec_inv,NBlengths);
    Monomials_simple[:] = 0;

    file = open("Invariant_generation/NBody_$NBody""_deg_$Deg""_invariants.jl")
    line = readlines(file)

    for i=1:length(line)
        for j=1:NBlengths
            if contains(line[i], "x[$j]")
                if contains(line[i], "x[$j]^")
                    for k=1:Deg
                        if contains(line[i], "x[$j]^$k")
                            Monomials[i,j] = k
                            Monomials_simple[i,j] = 1
                        end
                    end
                else
                    Monomials[i,j] = 1
                    Monomials_simple[i,j] = 1
                end
            end
        end
    end
    return NB_sec_inv,Monomials,Monomials_simple
end


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


function generate_file_1_perm(Perms,Permsref,filename,pref)
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
end




function generate_all_irr_sec(NBody,Deg)
   (NB_sec_inv,Monomials,Monomials_simple) = generate_monomials_irr_sec(NBody,Deg)
   filename = "NBody_$NBody"*"irr_sec_deg_$Deg.jl";
   open(filename, "w") do f
      write(f, "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n")
   end
   for j=1:NB_sec_inv
      pref = "PV$j";
      Perms = SVector(Monomials[j,:]...)
      Permsref = SVector(Monomials_simple[j,:]...)
      generate_file_1_perm(Perms,Permsref,filename,pref)
   end
   return 0
end


NBody = 4;
Deg = 10;
generate_all_irr_sec(NBody,Deg)
