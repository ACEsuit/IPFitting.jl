using Combinatorics, StaticArrays



# miscallaneous
Base.ntuple(n::Integer, ::Val{M}) where {M} = ntuple(_->n, M)


# --------------- helpful functions to connect ---------------
#                 permutations of corners and edges


# 4-body = 4-simplex has 4 corners (atom positions) and 6 edges  rᵢⱼ
#
# edge index : A[1]  A[2]  A[3]  A[4]  A[5]  A[6]
# edge length: r12   r13   r14   r23   r24   r34
# (where rij = |xᵢ - xⱼ| with xᵢ the corner positions)

const b4_e_inds = [0 1 2 3
                   1 0 4 5
                   2 4 0 6
                   3 5 6 0]


# 5-body = 5-simplex has 5 corners (atom positions) and 10 edges  rᵢⱼ
#
# edge index : A[1]  A[2]  A[3]  A[4]  A[5]  A[6] A[7]  A[8]  A[9]  A[10]
# edge length: r12   r13   r14   r15   r23   r24  r25   r34   r35   r45
# (where rij = |xᵢ - xⱼ| with xᵢ the corner positions)

const b5_e_inds = [0 1 2 3 4
                   1 0 5 6 7
                   2 5 0 8 9
                   3 6 8 0 10
                   4 7 9 10 0]


# preparation for true n-body terms only
 #for 3-body
 const b3_e_proj = [1 1 0
                    1 0 1
                    0 1 1]
#for 4-body
 const b4_e_proj = [1 1 1 0 0 0
                    1 0 0 1 1 0
                    0 1 0 1 0 1
                    0 0 1 0 1 1]
#for 5-body
const b5_e_proj = [1 1 1 1 0 0 0 0 0 0
                  1 0 0 0 1 1 1 0 0 0
                  0 1 0 0 1 0 0 1 1 0
                  0 0 1 0 0 1 0 1 0 1
                  0 0 0 1 0 0 1 0 1 1]

const πb3 = collect(permutations(1:3))

"""
convert a permutation of simplex corners into a permutation of
simplex edges (for 4 body)
"""
S4_to_S6(π::Vector{Int}, b4_e_inds=NBodyIPs.b4_e_inds) = Int[
   b4_e_inds[π[1], π[2]], b4_e_inds[π[1], π[3]], b4_e_inds[π[1], π[4]],
   b4_e_inds[π[2], π[3]], b4_e_inds[π[2], π[4]], b4_e_inds[π[3], π[4]] ]

"""
convert a permutation of simplex corners into a permutation of
simplex edges (for 5 body)
"""
S5_to_S10(π::Vector{Int}, b5_e_inds=NBodyIPs.b5_e_inds) = Int[
   b5_e_inds[π[1], π[2]], b5_e_inds[π[1], π[3]], b5_e_inds[π[1], π[4]],
   b5_e_inds[π[1], π[5]], b5_e_inds[π[2], π[3]], b5_e_inds[π[2], π[4]],
   b5_e_inds[π[2], π[5]], b5_e_inds[π[3], π[4]], b5_e_inds[π[3], π[5]],
   b5_e_inds[π[4], π[5]] ]


"""
generate all permutations of A that correspond to permutations of corners,
then keep only the unique ones so that we don't double-count.
"""
simplex_permutations(::Val{5}, A) =
   (  [ A[S5_to_S10(πX)]
              for πX in permutations(1:5) ]  )

simplex_permutations(::Val{4}, A) =
   (  [ A[S4_to_S6(πX)]
              for πX in permutations(1:4) ]  )

simplex_permutations(::Val{3}, A) = ( [ A[π] for π in πb3 ] )

simplex_permutations(::Val{2}, A) = [ A ]


simplex_permutations(x::SVector{3}) =
   simplex_permutations(Val(3), x)

simplex_permutations(x::SVector{6}) =
   simplex_permutations(Val(4), x)

simplex_permutations(x::SVector{10}) =
   simplex_permutations(Val(5), x)
