struct Vertex
    tau::Float64
    spin::Float64
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


    function QMC_variables(count, vertices, averaged_sign, P, N, Dr, Dc, S, Scount, mkink, norbs)
        R_temp = zeros(mkink)
        L_temp = zeros(mkink)
        Mklm_temp = zeros(mkink, mkink, norbs)
        ev_temp = zeros(Float64, norbs)
        g0l_temp = zeros(Float64, mkink, norbs)
        ntemp_temp = zeros(Float64, mkink, mkink)
        Mkl_temp = zeros(norbs)
        return new(count, vertices, averaged_sign, P, N, Dr, Dc, S, Scount, R_temp, L_temp, Mklm_temp, ev_temp, g0l_temp, ntemp_temp, Mkl_temp)
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

end

struct QMC
    v::QMC_variables
    p::QMC_parameters
end