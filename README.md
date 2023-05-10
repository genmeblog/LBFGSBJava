# LBFGSBJava

L-BFGS-B box constrained optimizer 
* Most of the code ported from [LBFGS++](https://lbfgspp.statr.me/) code.
* MoreThuente line search ported from [R source](https://github.com/SurajGupta/r-source/blob/master/src/appl/lbfgsb.c#L2976) and [Julia LineSearches.jl](https://github.com/JuliaNLSolvers/LineSearches.jl/blob/master/src/morethuente.jl)
* MoreThuente can be used with strong or weak (default) Wolfe condition 
* LewisOverton (experimental, don't use) line search ported from [LBFGS-Lite](https://github.com/ZJU-FAST-Lab/LBFGS-Lite/blob/master/include/lbfgs.hpp)

LBFGSBJava is implemented without any dependencies.

See `org.generateme.lbfgsb.examples` for usage cases.

## Target function

Each target function should implement `IGradFunction` interface.

* default implementation of `gradient` uses finite differences method
* implementing only `evaluate(x)` will use finite differences
* there is a possibility to calculate function value and gradient in one call, implement `gradient(x,grad)` and implement `in_place_gradient` to return `true`. See `Rosenbrock` in examples.

## Licence

Copyright (c) 2023 GenerateMe

The MIT Licence
