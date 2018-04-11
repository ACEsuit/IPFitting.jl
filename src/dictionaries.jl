export dict


dict(::Val{:poly}, n) =
    ["r^$i" for i = 1:n], "r"

dict(v::Val{:poly}, n, rcut) = dict(v, n)

dict(::Val{:poly1}, n, rcut) =
    ["x^$i * (x^(-1) - $(1/rcut) + $(1/rcut^2) * (x-$rcut))" for i = 0:n-1], "x"

dict(::Val{:poly2}, n, rcut) =
    ["(x*$(1/rcut)-1.0)^$(2+i)" for i = 0:n-1], "x"

dict(::Val{:inv1}, n, rcut) =
    ["x^$(-i) - $(rcut^(-i)) + $(i * rcut^(-i-1)) * (x - $rcut)" for i = 1:n], "x"

dict(::Val{:inv2}, n, rcut) =
    ["x^$(-i) * (x*$(1/rcut)-1.0)^2" for i = 0:n-1], "x"

dict(::Val{:exp1}, n, rcut) =
    ["exp(-$i * x) - $(exp(-i*rcut)) + $(i * exp(-i*rcut)) * (x-$rcut)" for i = 1:n], "x"

dict(sym, args...) = dict(Val(sym), args...)
