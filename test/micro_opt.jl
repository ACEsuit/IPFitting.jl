using StaticArrays, BenchmarkTools, Unrolled

const i1 = @SVector [1,1,1,1,1,1,2,2,3,4,3,4,2,2,3,4,3,4,2,2,3,4,3,4,5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,5,5,6,7,6,7,8,9,8,9,10,10,]
const i2 = @SVector [3,4,2,2,4,3,3,4,2,2,4,3,1,1,1,1,1,1,4,3,4,3,2,2,2,2,3,4,3,4,1,1,1,1,1,1,3,4,2,2,4,3,3,4,2,2,4,3,7,6,7,6,5,5,6,7,5,5,7,6,]
const i3 = @SVector [5,5,6,7,6,7,5,5,6,7,6,7,8,9,8,9,10,10,8,9,8,9,10,10,6,7,5,5,7,6,8,9,8,9,10,10,5,5,6,7,6,7,9,8,10,10,8,9,8,9,8,9,10,10,9,8,10,10,8,9,]
x = @SVector rand(10)

tdot_A(x) = (dot(x[i1] .* x[i2], x[i3]))
tdot_B(x) = sum(x[i1] .* x[i2] .* x[i3])
# tdot_C(x) = sum(x[i1[n]] .* x[i2[n]] .* x[i3[n]] for n = 1:60)
function tdot_D(x)
    y = zero(eltype(x))
    for n = 1:length(i1)
        @inbounds a = x[i1[n]] * x[i2[n]]
        @inbounds y = muladd(a, x[i3[n]], y)
    end
    return y
end

@btime tdot_A($x)
@btime tdot_B($x)
# @btime tdot_C($x)
@btime tdot_D($x)

# @code_llvm(tdot_A(x))
@code_llvm(tdot_D(x))

@btime ForwardDiff.gradient($tdot_A, $x)
@btime ForwardDiff.gradient($tdot_D, $x)

const j1 = @SVector [1,1,1,2,2,2,3,4,3,4,4,3,4,5,6,7,5,6,7,7,8,9,8,9,10,10,9,5,5,6,7,6,7,8,9,8,9,10,10,]
const j2 = @SVector [3,4,2,1,3,4,2,2,4,3,1,1,1,2,3,4,1,1,1,1,3,4,2,2,4 ,3 ,2,7,6,7,6,5,5,6,7,5,5,7 ,6 ,]
A12_ = rand([0,1], (10,10))
const A12 = SMatrix{10,10}(A12_...)

ddot_A(x) = (dot(x[j1], x[j2]))
ddot_B(x) = (dot(x, A12 * x) )
@inline @unroll function ddot_C(x, j1,j2)
   y = zero(eltype(x))
   @unroll for n = 1:length(x)
     @inbounds y = muladd(x[j1[n]], x[j2[n]], y)
   end
   return y
end

x = @SVector rand(10)

@btime ddot_A($x)
@btime ddot_B($x)
@btime ddot_C($x, $j1, $j2)
