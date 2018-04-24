using NBodyIPs
using JuLIP, Base.Test, StaticArrays
using BenchmarkTools

@testset "NBodyIPs" begin
   @testset "Invariants" begin include("test_invariants.jl") end
end

# TODO  - write test set for fitting

using StaticArrays, BenchmarkTools

@inline _m1d(α::Number, x::T) where {T <: Number} =
      α == 0 ? zero(T) : (@fastmath α * x^(α-1))

@inline _m1(α::Number, x::T) where {T <: Number} =
      (@fastmath x^α)

function _monomial3_d(α::NTuple{M, TI},
                             x::SVector{NI, T}) where {M, TI, NI, T}
   f1 = @fastmath x[1]^α[1]
   f2 = @fastmath x[2]^α[2]
   f3 = @fastmath x[3]^α[3]
   m = @fastmath f1 * f2 * f3
   m_d = SVector{4, T}( _m1d(α[1], x[1]) * f2 * f3,
                        _m1d(α[2], x[2]) * f1 * f3,
                        _m1d(α[3], x[3]) * f1 * f2, 0.0 )
   return m, m_d
end


   # SVector{4, T}( df[1] *  f[2] *  f[3],
   #                       f[1] * df[2] *  f[3],
   #                       f[1] *  f[2] * df[3], 0.0 )

function _pd(f::SVector{K,T}, df::SVector{K,T}) where {K, T}
   for n = 2:K
      f = [shift(f); SVector{1,T}(f[1])]
      df = df .* f
   end
   return df
end

function _monomial_d(α::NTuple{K, TI},
                     x::SVector{NI, T}, ::Val{3}) where {K, TI, NI, T}
   # f = @SVector [_m1(α[n], x[n]) for n = 1:M]
   # df = @SVector [_m1d(α[n], x[n]) for n = 1:M]
   f = SVector{3,T}(_m1(α[1], x[1]), _m1(α[2], x[2]), _m1(α[3], x[3]))
   df = SVector{3,T}(_m1d(α[1], x[1]), _m1d(α[2], x[2]), _m1d(α[3], x[3]))
   m = prod(f)
   m_d = _pd(f, df)
   return m, m_d
end

function monomial(α::NTuple, x::SVector{NN,T}, ::Val{M}) where {M, NN, T}
   m = one(T)
   for i = 1:M
      m *= @fastmath x[i]^α[i]
   end
   return m
end


@generated function monomial_gen(α::NTuple{K,TI},
                                   x::SVector) where {K, TI}
   ex = "@fastmath "
   for i = 1:K-1
      ex *= "x[$i]^α[$i] * "
   end
   ex = ex[1:end-3]
   quote
      $(parse(ex))
   end
end

@generated function monomial_d_gen(α::NTuple{K,TI},
                                   x::SVector{NX,T}) where {K, TI, NX,T}

   # evaluate the scalar monomials
   ex_f = "f = SVector{$(K-1), $T}("
   ex_df = "df = SVector{$(K-1), $T}("
   for i = 1:K-1
      ex_f  *= " _m1(α[$i], x[$i]), "
      ex_df *= "_m1d(α[$i], x[$i]), "
   end
   ex_f  =  ex_f[1:end-2] * ")"
   ex_df = ex_df[1:end-2] * ")"
   # @show ex_f
   # @show ex_df

   # evaluate the multi-variate monomial
   ex_m = "m = "
   for i = 1:K-1
      ex_m *= "f[$i] * "
   end
   ex_m = ex_m[1:end-3]
   # @show ex_m

   # evaluate the derivative
   ex_dm = "m_d = SVector{$(K-1), $T}("
   for i = 1:K-1
      dmj = ""
      for j = 1:K-1
         if i == j
            dmj *= " df[$j] *"
         else
            dmj *= "  f[$j] *"
         end
      end
      ex_dm *= dmj[1:end-2] * ", "
   end
   ex_dm = ex_dm[1:end-2] * ")"
   # @show ex_dm

   quote
      $(parse(ex_f))
      $(parse(ex_df))
      $(parse(ex_m))
      $(parse(ex_dm))
      return m, m_d
   end
end

      # $ex_m;
      # $ex_dm;

monomial_d_gen(α, x)

α = (3,5,1,0)
x = @SVector rand(4)
vN = Val(3)

@btime monomial($α, $x, $vN)
@btime monomial_gen($α, $x)

@btime _monomial3_d($α, $x)
@btime _monomial_d($α, $x, $vN)
@btime monomial_d_gen($α, $x)
