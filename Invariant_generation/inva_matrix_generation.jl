
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
