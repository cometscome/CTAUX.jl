
# CTAUX.jl: Continuous-time auxiliary-field Quantum Monte Carlo method for Anderson impurity model with Bethe-lattice-bath electrons

See, E. Gull et al., EPL 82, 57003 (2008)

[![Build Status](https://github.com/cometscome/CTAUX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/cometscome/CTAUX.jl/actions/workflows/CI.yml?query=branch%3Amain)

# install

```
add https://github.com/cometscome/CTAUX.jl 
```

# example

```julia
using CTAUX
function test()
    beta = 10.0
    U = 2.0
    mu = U / 2
    K = 1.0
    V = 1.0 #Strength of the hybridization

    system = QMC(U, K, beta, V, mu)
    #display(system)
    nthermal = 1000


    mqs = 1000000
    run_ctaux!(system, nthermal, measurements=false)
    run_ctaux!(system, mqs, measurements=true)

    println("Calculating Green's function...")
    Gτ = calc_green(system)
    tau = get_taumesh(system)
end
test()
```