struct Vertex
    tau::Float64
    spin::Int64#Float64
end

mutable struct Vertices
    num::Int64
    values::Array{Vertex,1}
    indices::Array{Int64,1}
    currentk::Int64
end

struct Greenfunction
    spl::Dierckx.Spline1D
end

function (g::Greenfunction)(x)
    return g.spl(x)
end


function get_verticesinfo(v::Vertices)
    taus = Tuple{Float64,Int64}[]
    for i = 1:v.currentk
        index = v.indices[i]
        spin = v.values[index].spin
        tau = v.values[index].tau
        push!(taus, (tau, spin))
    end
    return taus
end



mutable struct QMC_variables
    count::Int64
    vertices::Vertices
    averaged_sign::Float64
    P::Array{Int64,1} #Distribution of the orders
    N::Array{Float64,3}
    Dr::Array{Float64,2}
    Dc::Array{Float64,2}
    S::Array{Float64,2}
    Scount::Array{Float64,2}
    R_temp::Vector{Float64}
    L_temp::Vector{Float64}
    Mklm_temp::Array{Float64,3}
    ev_temp::Vector{Float64}
    g0l_temp::Matrix{Float64}
    ntemp_temp::Matrix{Float64}
    Mkl_temp::Vector{Float64}
    current_logW::Float64


    function QMC_variables(count, vertices, averaged_sign, P, N, Dr, Dc, S, Scount, mkink, norbs)
        R_temp = zeros(mkink)
        L_temp = zeros(mkink)
        Mklm_temp = zeros(mkink, mkink, norbs)
        ev_temp = zeros(Float64, norbs)
        g0l_temp = zeros(Float64, mkink, norbs)
        ntemp_temp = zeros(Float64, mkink, mkink)
        Mkl_temp = zeros(norbs)
        current_logW = 0
        return new(count, vertices, averaged_sign, P, N, Dr, Dc, S, Scount,
            R_temp, L_temp, Mklm_temp, ev_temp,
            g0l_temp, ntemp_temp, Mkl_temp,
            current_logW)
    end
end

struct QMC_parameters
    ntime::Int64
    mfreq::Int64
    mu::Float64
    K::Float64
    U::Float64
    beta::Float64
    mesh_tau::Array{Float64,1}
    mesh_omega::Array{Float64,1}
    L::Int64
    norb::Int64
    expgamma0::Tuple{Float64,Float64}
    delta::Array{ComplexF64,3}
    gtau0::Array{Greenfunction,1}


    function QMC_parameters(beta, U, K, mu, V; ntime=1024, mfreq=1024, L=1, norb=2)

        gamma0 = acosh(1.0 + (beta * U) / (2.0 * K))
        expgamma0 = (exp(gamma0), exp(-gamma0))
        mesh_tau = calc_linear_mesh(0, beta, ntime)
        mesh_omega = calc_linear_mesh(π / beta, (π / beta) * (2mfreq - 1.0), mfreq)
        Ef = -mu

        delta = make_hybridization(mesh_omega, V, norb)

        Gf0_omega = make_Gf0_omega(mesh_omega, norb, Ef, delta, U)

        Gf0_tau = make_Gf0_tau(Gf0_omega, mesh_omega, ntime, beta)
        Gf0_tau = -Gf0_tau

        gtau0 = set_greenfunctions(Gf0_tau, mesh_tau)

        return new(ntime, mfreq, mu, K, U, beta, mesh_tau, mesh_omega, L, norb, expgamma0, delta, gtau0)
    end
end

