using Combinatorics, StaticArrays, NBodyIPs

# include("misc.jl")



# TODO: code this function in a nicer way
function generate_monomials_irr_sec(NBody,Deg)
    NBlengths = Int(NBody*(NBody-1)/2);
    NB_sec_inv = countlines("Invariant_generation/NBody_$NBody""_deg_$Deg""_invariants.jl");

    Monomials = [];
    Monomials_simple = [];

    file = open("Invariant_generation/NBody_$NBody""_deg_$Deg""_invariants.jl")
    line = readlines(file)

    for i=1:length(line)
        Mono_temp = 0*Vector{Int64}(NBlengths)
        Mono_sim_temp = 0*Vector{Int64}(NBlengths)
        for j=1:NBlengths
            if contains(line[i], "x[$j]")
                if contains(line[i], "x[$j]^")
                    for k=1:Deg
                        if contains(line[i], "x[$j]^$k")
                            Mono_temp[j] = k
                            Mono_sim_temp[j] = 1
                        end
                    end
                else
                    Mono_temp[j] = 1
                    Mono_sim_temp[j] = 1
                end
            end
        end
        push!(Monomials,Mono_temp);
        push!(Monomials_simple,Mono_sim_temp);
    end
    return NB_sec_inv,Monomials,Monomials_simple
end

generate_monomials_irr_sec(4,10)

function printmatrix(filename,A)
    M = size(A,1);
    N = size(A,2);
    open(filename, "a") do f
       write(f, "[")
       for m=1:(M-1)
           for n=1:N
               element = A[m,n]
                write(f, "$element")
                write(f, " ")
            end
            write(f, "; ")
        end
            for n=1:N
            element = A[M,n]
             write(f, "$element")
             write(f, " ")
         end
        write(f, "]\n")
    end
end

function _generate_matrix_Deg_2(NBody,Perms)
    NBlengths = Int(NBody*(NBody-1)/2);
    Ind = perm_2_indice(unique(simplex_permutations(Perms)))
    A = zeros(Int64, NBlengths,NBlengths);
    for i=1:size(Ind, 1)
        A[Ind[i,1],Ind[i,2]] = 1;
    end
    return A
end

A= _generate_matrix_Deg_2(4,SVector(2,1,0,0,0,0))



function generate_all_irr_sec(NBody,Deg)
   (NB_sec_inv,Monomials,Monomials_simple) = generate_monomials_irr_sec(NBody,Deg)
   filename = "Invariant_generation/NBody_$NBody"*"irr_sec_deg_$Deg.jl";
   open(filename, "w") do f
      write(f, "# Irreducible secondaries for NBody=$NBody"*"and deg=$Deg \n")
   end
   for j=1:NB_sec_inv
      pref = "PV$j";
      Perms = SVector(Monomials[j]...)
      Permsref = SVector(Monomials_simple[j]...)
      perm_deg = sum(Permsref)
      nb_diff_coef = length(unique(Perms))-1
      Perm_inv = sort(Perms,rev=true)
      # TODO: keep track of matrices

      # different case with respect to the number of variables in the monomial
      if (perm_deg==1)
          # only 1 variable -> sum of x to some power
          open(filename, "a") do f
             exponent = Perm_inv[1];
             write(f, pref)
             write(f, " = ")
             write(f, "sum(x$exponent)\n")
          end

      elseif (perm_deg==2)
          A = _generate_matrix_Deg_2(NBody,Perms)
          exp1 = Perm_inv[1];
          exp2 = Perm_inv[2];
          open(filename, "a") do f
             write(f, "A$j=")
         end
             printmatrix(filename,A)
         open(filename, "a") do f
             write(f,"\n")
             write(f, pref)
             write(f, " = ")
             write(f, "x$exp1","'*","A$j","*"," x$exp2")
             write(f, "\n\n")
          end

      elseif (perm_deg==3)
          for s=1:perm_deg
          exp = Perm_inv[s]
          Ind = perm_2_indice(unique(simplex_permutations(Perms)))
          exp1 = Perm_inv[1];
          exp2 = Perm_inv[2];
          exp3 = Perm_inv[3];
          open(filename, "a") do f
              # First vector
              write(f, pref, "_1 = ")
              write(f, "@SVector [")
              for j = 1:size(Ind,1)-1
                  k = Ind[j,1]
                  write(f, "x$exp1[$k],")
              end
              write(f, "] \n")
              # Second vector
              write(f, pref, "_2 = ")
              write(f, "@SVector [")
              for j = 1:size(Ind,1)-1
                  k = Ind[j,2]
                  write(f, "x$exp2[$k],")
              end
              write(f, "] \n")
              # Third vector
              write(f, pref, "_3 = ")
              write(f, "@SVector [")
              for j = 1:size(Ind,1)-1
                  k = Ind[j,3]
                  write(f, "x$exp3[$k],")
              end
              write(f, "] \n")
              write(f, pref, " = ")
              write(f, "dot(", pref, "_1.*", pref, "_2,", pref, "_3)\n \n")
         end

      elseif (perm_deg==4)

      elseif (perm_deg==5)

      else
          error("not implemented (perm_deg too large)")
      end
   end
   return 0
end


generate_all_irr_sec(5,6)

function perm_2_indice(Perms)
   L = length(Perms);
   deg = length(find(Perms[1]));
   M = length(Perms[1]);

   Indices = Array{Int64, 2}(L,deg);
   # zeros(L,deg);

      for j=1:L
         Ind_temp = sortperm(Perms[j],rev=true)
         Indices[j,:] = Ind_temp[1:deg]
      end
   return Indices
end

perm_2_indice([SVector(2,1,0,0,0,0),SVector(1,0,0,0,0,2)])
2;

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
   elseif (Decreasing_coef[1]==2)&(Decreasing_coef[2]==1)&(perm_deg==2)
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
   elseif (Decreasing_coef[1]==4)&(Decreasing_coef[2]==1)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[1],"_4 =",Names[1],"_3.*",Names[1],"\n")
         write(f,Names[2],"_4 =",Names[2],"_3.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_4,"*Names[2]*")+dot("*Names[1]","Names[2]"_4)\n")
      end
   elseif (Decreasing_coef[1]==5)&(Decreasing_coef[2]==1)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[1],"_4 =",Names[1],"_3.*",Names[1],"\n")
         write(f,Names[2],"_4 =",Names[2],"_3.*",Names[2],"\n")
         write(f,Names[1],"_5 =",Names[1],"_4.*",Names[1],"\n")
         write(f,Names[2],"_5 =",Names[2],"_4.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_5,"*Names[2]*")+dot("*Names[1]","Names[2]"_5)\n")
      end
   elseif (Decreasing_coef[1]==6)&(Decreasing_coef[2]==1)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[1],"_4 =",Names[1],"_3.*",Names[1],"\n")
         write(f,Names[2],"_4 =",Names[2],"_3.*",Names[2],"\n")
         write(f,Names[1],"_5 =",Names[1],"_4.*",Names[1],"\n")
         write(f,Names[2],"_5 =",Names[2],"_4.*",Names[2],"\n")
         write(f,Names[1],"_6 =",Names[1],"_5.*",Names[1],"\n")
         write(f,Names[2],"_6 =",Names[2],"_5.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_6,"*Names[2]*")+dot("*Names[1]","Names[2]"_6)\n")
      end
   elseif (Decreasing_coef[1]==2)&(Decreasing_coef[2]==2)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_2,"*Names[2]*"_2)\n")
      end
   elseif (Decreasing_coef[1]==2)&(Decreasing_coef[2]==1)&(Decreasing_coef[3]==1)&(perm_deg==3)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_2.*"Names[2]","*Names[3]")+dot("*Names[1]*".*"Names[2]"_2,"*Names[3]")+dot("*Names[1]*".*"Names[2]","*Names[3]"_2)\n")
      end
   elseif (Decreasing_coef[1]==3)&(Decreasing_coef[2]==2)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_2,"Names[2]"_3)+dot("*Names[1]*"_3,"Names[2]"_2)\n")
      end
   elseif (Decreasing_coef[1]==3)&(Decreasing_coef[2]==1)&(Decreasing_coef[3]==1)&(perm_deg==3)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[3],"_3 =",Names[3],"_2.*",Names[3],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_3.*"Names[2]","*Names[3]")+dot("*Names[1]*".*"Names[2]"_3,"*Names[3]")+dot("*Names[1]*".*"Names[2]","*Names[3]"_3)\n")
      end
   elseif (Decreasing_coef[1]==2)&(Decreasing_coef[2]==2)&(Decreasing_coef[3]==1)&(perm_deg==3)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_2.*"Names[2]"_2,"*Names[3]")+dot("*Names[1]*"_2.*"Names[2]","*Names[3]"_2)+dot("*Names[1]*".*"Names[2]"_2,"*Names[3]"_2)\n")
      end
   elseif (Decreasing_coef[1]==4)&(Decreasing_coef[2]==2)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[1],"_4 =",Names[1],"_3.*",Names[1],"\n")
         write(f,Names[2],"_4 =",Names[2],"_3.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_4,"*Names[2]*"_2)+dot("*Names[1]"_2,"Names[2]"_4)\n")
      end
   elseif (Decreasing_coef[1]==3)&(Decreasing_coef[2]==3)&(perm_deg==2)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_3,"*Names[2]*"_3)\n")
      end
   elseif (Decreasing_coef[1]==4)&(Decreasing_coef[2]==1)&(Decreasing_coef[3]==1)&(perm_deg==3)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[3],"_3 =",Names[3],"_2.*",Names[3],"\n")
         write(f,Names[1],"_4 =",Names[1],"_3.*",Names[1],"\n")
         write(f,Names[2],"_4 =",Names[2],"_3.*",Names[2],"\n")
         write(f,Names[3],"_4 =",Names[3],"_3.*",Names[3],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_4.*"Names[2]","*Names[3]")+dot("*Names[1]*".*"Names[2]"_4,"*Names[3]")+dot("*Names[1]*".*"Names[2]","*Names[3]"_4)\n")
      end
   elseif (Decreasing_coef[1]==3)&(Decreasing_coef[2]==2)&(Decreasing_coef[3]==1)&(perm_deg==3)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[3],"_3 =",Names[3],"_2.*",Names[3],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*".*"Names[2]"_2,"*Names[3]"_3)+dot("*Names[1]*".*"Names[2]"_3,"*Names[3]"_2)+dot("*Names[1]*"_2.*"Names[2]","*Names[3]"_3) +dot("*Names[1]*"_2.*"*Names[2]*"_3,"*Names[3]*")+dot("*Names[1]*"_3.*"Names[2]","*Names[3]"_2)+dot("*Names[1]*"_3.*"Names[2]"_2,"*Names[3]*")\n")
      end
   elseif (Decreasing_coef[1]==2)&(Decreasing_coef[2]==2)&(Decreasing_coef[3]==2)&(perm_deg==3)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_2.*"*Names[2]*"_2,"*Names[3]*"_2)\n")
      end
   elseif (Decreasing_coef[1]==3)&(Decreasing_coef[2]==1)&(Decreasing_coef[3]==1)&(Decreasing_coef[4]==1)&(perm_deg==4)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f,Names[4],"_2 =",Names[4],".*",Names[4],"\n")
         write(f,Names[1],"_3 =",Names[1],"_2.*",Names[1],"\n")
         write(f,Names[2],"_3 =",Names[2],"_2.*",Names[2],"\n")
         write(f,Names[3],"_3 =",Names[3],"_2.*",Names[3],"\n")
         write(f,Names[4],"_3 =",Names[4],"_2.*",Names[4],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*"_3.*"*Names[2]*","*Names[3]*".*"*Names[4]*")
         +dot("*Names[1]*".*"*Names[2]*"_3,"*Names[3]*".*"*Names[4]*") +dot("*Names[1]*".*"*Names[2]*","*Names[3]*"_3.*"*Names[4]*") +dot("*Names[1]*".*"*Names[2]*","*Names[3]*".*"*Names[4]*"_3)\n")
      end
   elseif (Decreasing_coef[1]==2)&(Decreasing_coef[2]==2)&(Decreasing_coef[3]==1)&(Decreasing_coef[4]==1)&(perm_deg==4)
      open(filename, "a") do f
         write(f,Names[1],"_2 =",Names[1],".*",Names[1],"\n")
         write(f,Names[2],"_2 =",Names[2],".*",Names[2],"\n")
         write(f,Names[3],"_2 =",Names[3],".*",Names[3],"\n")
         write(f,Names[4],"_2 =",Names[4],".*",Names[4],"\n")
         write(f, pref)
         write(f, " = ")
         write(f, "dot("*Names[1]*".*"*Names[2]*","*Names[3]*"_2.*"*Names[4]*"_2)    +dot("*Names[1]*".*"*Names[2]*"_2,"*Names[3]*"_2.*"*Names[4]*")
         +dot("*Names[1]*".*"*Names[2]*"_2,"*Names[3]*".*"*Names[4]*"_2)         +dot("*Names[1]*"_2.*"*Names[2]*","*Names[3]*".*"*Names[4]*"_2) +dot("*Names[1]*"_2.*"*Names[2]*","*Names[3]*"_2.*"*Names[4]*") +dot("*Names[1]*"_2.*"*Names[2]*"_2,"*Names[3]*".*"*Names[4]*") \n")
      end
   else
      error("not implemented yet (perm_svector_generator)")
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
