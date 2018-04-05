# NBodyIPs.jl

<!--
[![Build Status](https://travis-ci.org/cortner/NBodyIPs.jl.svg?branch=master)](https://travis-ci.org/cortner/NBodyIPs.jl)

[![Coverage Status](https://coveralls.io/repos/cortner/NBodyIPs.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/cortner/NBodyIPs.jl?branch=master)

[![codecov.io](http://codecov.io/github/cortner/NBodyIPs.jl/coverage.svg?branch=master)](http://codecov.io/github/cortner/NBodyIPs.jl?branch=master)
-->


## Examples

```
PermPolys(3, ["x^0","x^1","x^2", "x^3"], "x")
```
generates the expressions and associated functions
```
:(x[1]^0)
:(x[1]^0 * x[2]^0 * x[3]^1 + x[1]^0 * x[2]^1 * x[3]^0 + x[1]^1 * x[2]^0 * x[3]^0)
:(x[1]^0 * x[2]^0 * x[3]^2 + x[1]^0 * x[2]^2 * x[3]^0 + x[1]^2 * x[2]^0 * x[3]^0)
:(x[1]^0 * x[2]^1 * x[3]^1 + x[1]^1 * x[2]^0 * x[3]^1 + x[1]^1 * x[2]^1 * x[3]^0)
:(x[1]^0 * x[2]^0 * x[3]^3 + x[1]^0 * x[2]^3 * x[3]^0 + x[1]^3 * x[2]^0 * x[3]^0)
:(x[1]^0 * x[2]^1 * x[3]^2 + x[1]^0 * x[2]^2 * x[3]^1 + x[1]^1 * x[2]^0 * x[3]^2 + x[1]^1 * x[2]^2 * x[3]^0 + x[1]^2 * x[2] ^0 * x[3]^1 + x[1]^2 * x[2]^1 * x[3]^0)
:(x[1]^1 * x[2]^1 * x[3]^1)
```
