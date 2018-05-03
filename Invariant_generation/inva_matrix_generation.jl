using Combinatorics, StaticArrays

include("misc.jl")
include("perm_svector_generator.jl")

NBody = 5;
Deg = 6;
NBlengths = Int(NBody*(NBody-1)/2);

NB_sec_inv = countlines("Invariant_generation/NBody_$NBody""_deg_$Deg""_invariants.jl");

Monomials = Array{Int64, 2}(NB_sec_inv,NBlengths);
Monomials[:] = 0;

file = open("Invariant_generation/NBody_$NBody""_deg_$Deg""_invariants.jl")
line = readlines(file)


counter = 0;
for i=1:length(line)
    for j=1:NBlengths
        if contains(line[i], "x[$j]")
            if contains(line[i], "x[$j]^")
                for k=1:Deg
                    if contains(line[i], "x[$j]^$k")
                        Monomials[i,j] = k
                    end
                end
            else
                Monomials[i,j] = 1
            end
        end
    end
end
display(Monomials)

# open("Invariant_generation/NBody_4_deg_10_invariants.jl") do f
#     for i in enumerate(eachline(f))
#       println(i)
#       str = parse(Float64,i)
#       display(str)
#     end
# end
