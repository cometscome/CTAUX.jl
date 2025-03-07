using CTAUX
using Test
using Plots

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
    vertices = get_verticesinfo(system)
    display(vertices)

    println("Calculating Green's function...")
    Gτ = calc_green(system)
    tau = get_taumesh(system)

    orderdisp = system.v.P
    #display(orderdisp)

    #display(Gr)

    filename = "greenfromfortran.txt"
    fd = open(filename, "r")
    times = zeros(Float64, 1024)
    Gtau = zeros(Float64, 1024, 2)

    cnt = 0
    while !eof(fd)
        cnt += 1
        line = readline(fd)
        u = split(line)

        times[cnt] = parse(Float64, u[1])
        Gtau[cnt, 1] = parse(Float64, u[2])
        Gtau[cnt, 2] = parse(Float64, u[3])
    end
    close(fd)


    plot(tau[:], Gτ[:, 1], label="Julia", title="Imaginary-time Green's function")
    plot!(tau[:], Gtau[:, 1], label="Fortran")
    savefig("comparison.png")




    plot(0:29, orderdisp[1:30] / sum(orderdisp[:]))
    savefig("order.png")

end

@testset "CTAUX.jl" begin
    # Write your tests here.
    test()
end
