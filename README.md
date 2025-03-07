
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

## configuration format
We can output configuration files. 

```julia
function test()
    beta = 10.0
    U = 2.0
    mu = U / 2
    K = 1.0
    V = 1.0 #Strength of the hybridization

    system = QMC(U, K, beta, V, mu)
    #display(system)
    nthermal = 1000


    datafilename = "config.txt"

    mqs = 1000000
    run_ctaux!(system, nthermal, measurements=false)
    run_ctaux!(system, mqs, measurements=true, recordconfig=true,
        recordinterval=100, datafilename=datafilename)

    println("Calculating Green's function...")
    Gτ = calc_green(system)
    tau = get_taumesh(system)
end
test()
```

Data format is shown as follows. 
```
26.078212443584423	10	2.5033917166631756	1	3.4788601691684917	1	4.369199082405516	1	4.740772623176241	1	6.131248187371384	1	6.472244188566955	1	7.021683242783143	1	7.76620382318522	1	9.713463484912184	1	9.995111315676402	1		
```
The first double number is $log W$. The second integer is a number of spins on the imaginary time axis. Others are written as tau_1 s_1 tau_2 s_2...