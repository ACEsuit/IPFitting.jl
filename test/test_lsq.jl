# TODO: review and rewrite this

using ACE1, JuLIP, Test, IPFitting
using JuLIP.Potentials: evaluate_d
using IPFitting: Dat, LsqDB, observation, vec_obs
const Lsq = IPFitting.Lsq
using LinearAlgebra: qr, norm

# generate random data
function generate_data(species, L, rmax, N, calc)
   data = []
   data_dat = Dat[]
   for n = 1:N
      at = bulk(species; cubic=true, pbc=true) * L
      rattle!(at, rand() * rmax)
      E = energy(calc, at)
      F = forces(calc, at)
      push!(data, (at = at, E = E, F = F))
      push!(data_dat, Dat(at, "rand"; E = E, F = F))
   end
   return data, data_dat
end

r0 = rnn(:Cu)
calc = let r0=r0
   LennardJones(r0=r0) * C2Shift(2.5*r0)
end
train, train_dat = generate_data(:Cu, 3, 0.33*r0, 50, calc)
rcut2 = cutoff(calc)

basis = rpi_basis(species = :Cu,
      N = 2,                       # correlation order = body-order - 1
      maxdeg = 10,            # polynomial degree
      r0 = r0,                      # estimate for NN distance
      #D = SparsePSHDegree(; wL=1.3, csp=1.0),
      rin = 0.5*r0, rcut = rcut2,   # domain for radial basis (cf documentation) #5.5
      pin = 2);

## generate a LSQ system by hand

function lsq(train)
   nobs = sum( 1+3*length(t.at) for t in train )
   Ψ = zeros(nobs, length(basis))
   y = zeros(nobs)

   irow = 0
   for (at, E, F) in train
      # add energy to the lsq system
      irow += 1
      y[irow] = E / length(at)
      Ψ[irow, :] = energy(basis, at) * 1.345  / length(at)

      # add forces to the lsq system
      nf = 3*length(at)
      y[(irow+1):(irow+nf)] = mat(F)[:]
      Fb = forces(basis, at)
      for ib = 1:length(basis)
         Ψ[(irow+1):(irow+nf), ib] = 2.987 * mat(Fb[ib])[:]
      end
      irow += nf
   end

   return Ψ, y
end

Ψ, Y = lsq(train)

# generate a database so that we generate the row indices
db = LsqDB("", basis, train_dat)

## test that Ψ_man, F_man are equivalent to a system generated from a db
Ψ_man, Y_man = Ψ, Y
Ψ, Y = IPFitting.Lsq.get_lsq_system( db,
               Vref = OneBody(:Cu => 0.0),
               weights = Dict("default" => Dict("E" => 1.345, "F" => 2.987)),
            )

##
println(@test Y ≈ Y_man)
println(@test Ψ ≈ Ψ_man)
