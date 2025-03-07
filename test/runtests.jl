using CTAUX
using Test
using Plots

function test()
    beta = 10.0
    U = 2.0
    U = 5.0
    mu = U / 2
    K = 1.0
    V = 1.0 #Strength of the hybridization

    system = QMC(U, K, beta, V, mu)
    #display(system)
    nthermal = 10000

    datafilename = "config.txt"

    mqs = 1000000
    run_ctaux!(system, nthermal, measurements=false)
    run_ctaux!(system, mqs, measurements=true, recordconfig=true,
        recordinterval=100, datafilename=datafilename)
    #vertices = get_verticesinfo(system)
    #display(vertices)

    data = readlines(datafilename)
    num = countlines(datafilename)
    logw = zeros(num)
    for i = 1:num
        logw[i] = parse(Float64, split(data[i])[1])
    end
    histogram(logw, bins=60)
    savefig("histogram.png")



    println("Calculating Green's function...")
    Gτ = calc_green(system)
    tau = get_taumesh(system)

    orderdisp = system.v.P
    #display(orderdisp)

    #display(Gr)
    #=

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

    =#


    plot(tau[:], Gτ[:, 1], label="Julia", title="Imaginary-time Green's function")
    savefig("green.png")
    #plot!(tau[:], Gtau[:, 1], label="Fortran")
    #savefig("comparison.png")




    plot(0:29, orderdisp[1:30] / sum(orderdisp[:]))
    savefig("order.png")
    @test true
end

@testset "CTAUX.jl" begin
    # Write your tests here.
    test()
end
