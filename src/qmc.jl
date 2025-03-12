

struct QMC
    v::QMC_variables
    p::QMC_parameters

    function QMC(U, K, beta, V, mu; ntime=1024, mfreq=1024, mkink=1024, L=1, seednum=1234, norb=2)
        count = 0
        Random.seed!(seednum)

        vertices = Vertices(0, Vertex[], Int64[], 0)

        p = QMC_parameters(beta, U, K, mu, V; ntime, mfreq, L, norb)
        #p = QMC_parameters(ntime, mfreq, mu, K, U, beta, mesh_tau, mesh_omega, L, norb, expgamma0, delta, gtau0)

        averaged_sign = 0
        P = zeros(Int64, mkink)
        N = zeros(Float64, mkink, mkink, norb)
        Dr = zeros(Float64, mkink, norb)
        Dc = zeros(Float64, mkink, norb)
        #global S,Scount
        S = zeros(Float64, ntime, norb)
        Scount = zeros(Float64, ntime, norb)

        v = QMC_variables(count, vertices, averaged_sign, P, N, Dr, Dc, S, Scount, mkink, norb)


        return new(v, p)
    end
end
export QMC

function print_steps(istep, numsteps, averagek, averagesign, number_accept)
    if istep % (numsteps ÷ 10) == 0
        println("Number of QMC steps: ", istep, " of ", numsteps)
        println("Average order k: ", averagek / istep)
        println("Average sign :", averagesign / istep)
        println("Total", "\t", "accepted", "\t", "rate")
        println(istep, "\t", number_accept, "\t", number_accept / istep)
    end
end


function run_ctaux!(q::QMC, numsteps; measurements=true, recordconfig=false, recordinterval=10000, datafilename="config.txt")
    number_accept = 0
    averagesign = 0.0
    averagek = 0.0
    if recordconfig
        fp = open(datafilename, "w")
    end

    for istep = 1:numsteps
        pass, rsign = update!(q)
        averagesign += rsign
        averagek += get_currentk(q)
        number_accept += ifelse(pass, 1, 0)

        print_steps(istep, numsteps, averagek, averagesign, number_accept)

        if measurements
            calc_S!(q)
        end

        if recordconfig && (istep % recordinterval) == 1
            vertices = get_verticesinfo(q)
            logW = q.v.current_logW
            numk = length(vertices)
            print(fp, logW, "\t", numk, "\t")
            for i = 1:numk
                print(fp, vertices[i][1], "\t", vertices[i][2], "\t")
            end
            println(fp, "\t")
        end

    end

    if recordconfig
        close(fp)
    end

    if measurements
        ntime = q.p.ntime
        S = get_S(q)
        S ./= numsteps
        Scount = get_Scount(q)
        Scount .*= ntime / sum(Scount)

        for i = 1:ntime
            S[i, :] ./= Scount[i]
        end
    end
    println("Average sign = ", averagesign / numsteps)
    return
end
export run_ctaux!

function get_verticesinfo(q::QMC)
    return get_verticesinfo(q.v.vertices)
end
export get_verticesinfo

function update!(q::QMC)
    r = rand()
    #println(r)

    if r > 0.5
        pass, rsign = qmc_insert!(q)
    else
        #            print_vertices(q)
        pass, rsign = qmc_remove!(q)

    end
    #println("pass $pass")
    return pass, rsign
end

function calc_green(q::QMC)
    S = get_S(q)

    τmesh = q.p.mesh_tau
    β = get_beta(q)
    ntime = get_ntime(q)
    norbs = get_norb(q)
    Gτ0spl = q.p.gtau0
    Gr = calc_green(S, τmesh, β, ntime, norbs, Gτ0spl)
    return Gr
end
export calc_green



function calc_S!(q::QMC)
    currentk = get_currentk(q)
    norbs = get_norb(q)
    Mklm = get_Mklm_temp(q, currentk)
    #Mklm = zeros(Float64,currentk,currentk,norbs)
    ev = get_ev_temp(q) #zeros(Float64,norbs)
    g0l = get_g0l_temp(q, currentk)#zeros(Float64,currentk,norbs)
    β = get_beta(q)
    ntime = get_ntime(q)
    S = get_S(q)
    Scount = get_Scount(q)
    Mkl = get_Mkl_temp(q)
    #global γ

    for k in 1:currentk
        kk = get_indices(q)[k]
        vertex = get_timeordered_vertex(q, k)
        #            kk = indexcon[k]
        τk = vertex.tau #τcon[kk]
        sk = vertex.spin #spincon[kk]
        for sigma in 1:norbs
            ssign = (-1)^(sigma - 1)
            ev[sigma] = ifelse(ssign * sk == 1, get_expgamma0(q, 1), get_expgamma0(q, 2))
            #ev[sigma] = exp(γ*ssign*sk)
        end
        for l in 1:currentk
            ll = get_indices(q)[l] #indexcon[l]
            τl = get_vertices(q).values[ll].tau #τcon[ll]
            for sigma in 1:norbs
                Mklm[kk, ll, sigma] = (ev[sigma] - 1) * get_N(q)[kk, ll, sigma]
            end
        end
    end

    for l in 1:currentk
        ll = get_indices(q)[l]
        #            ll = indexcon[l]
        τl = get_vertices(q).values[ll].tau
        #            τl = τcon[ll]
        for sigma in 1:norbs
            g0l[ll, sigma] = calc_Gf0_tau(q, τl, sigma) #Gτ0spl[sigma](τl)
        end
    end

    #        global ntime,β
    #Mkl = zeros(Float64,norbs)
    Mkl .= 0.0
    ξ = β / ntime
    id = 3
    for k in 1:currentk
        kk = get_indices(q)[k]
        #kk = indexcon[k]
        τk = get_vertices(q).values[kk].tau
        #τk = τcon[kk]
        iτ = ceil(Int64, ntime * τk / β)
        Mkl[:] .= 0.0#zeros(Float64,norbs)
        for l in 1:currentk
            ll = get_indices(q)[l]
            #                ll = indexcon[l]
            for sigma in 1:norbs
                Mkl[sigma] += Mklm[kk, ll, sigma] * g0l[ll, sigma]
            end
        end

        for j in -id:id
            jj = iτ + j
            if 1 <= jj <= ntime
                τd = get_tau(q, jj)
                #τd = τmesh[jj]
                f = gauss(ξ, τd, τk)
                for sigma in 1:norbs
                    S[jj, sigma] += f * Mkl[sigma]
                end
                Scount[jj] += f
            end
        end

    end


end




function gauss(ξ, x, x0)
    f = (1 / (sqrt(2π) * ξ)) * exp(-(x - x0)^2 / (2 * ξ^2))
    return f
end

@inline function get_tau(q::QMC, itau::Int)
    return @inbounds q.p.mesh_tau[itau]
end

function get_taumesh(q::QMC)
    return q.p.mesh_tau
end
export get_taumesh

@inline function get_timeordered_vertex(q::QMC, i)
    return @inbounds get_values(q)[get_indices(q)[i]]
end

function get_tau_from_vertices(q::QMC, i)
    return get_timeordered_vertex(q, i).tau
end

function get_spin_from_vertices(q::QMC, i)
    return get_timeordered_vertex(q, i).tau
end

@inline function get_omega(q::QMC, iomega::Int)
    return @inbounds q.p.mesh_omega[iomega]
end

#=
@inline function get_expgamma0(q::QMC, upordown)
    i = ifelse(upordown == 1, 1, 2)
    return @inbounds q.p.expgamma0[i]
end
=#

@inline function get_U(q::QMC)
    return q.p.U
end

@inline function get_K(q::QMC)
    return q.p.K
end

@inline function get_L(q::QMC)
    return q.p.L
end

@inline function get_norb(q::QMC)
    return q.p.norb
end

@inline function get_beta(q::QMC)
    return q.p.beta
end

@inline function get_ntime(q::QMC)
    return q.p.ntime
end

@inline function get_mfreq(q::QMC)
    return q.p.mfreq
end

@inline function get_count(q::QMC)
    return q.v.count
end

@inline function nextcount!(q::QMC)
    q.v.count += 1
    return
end

@inline function get_N(q::QMC)
    return q.v.N
end

@inline function get_S(q::QMC)
    return q.v.S
end

@inline function get_Scount(q::QMC)
    return q.v.Scount
end


@inline function get_Dr(q::QMC)
    return q.v.Dr
end

@inline function get_Dr(q::QMC, k, iorb)
    return view(q.v.Dr, 1:k, iorb)
end

@inline function get_Dc(q::QMC)
    return q.v.Dc
end

@inline function get_Dc(q::QMC, k, iorb)
    return view(q.v.Dc, 1:k, iorb)
end



@inline function get_vertices(q::QMC)
    return q.v.vertices
end

@inline function get_currentk(q::QMC)
    return get_vertices(q).currentk
end


@inline function currentk_add1!(q::QMC)
    get_vertices(q).currentk += 1
    return
end

@inline function currentk_remove1!(q::QMC)
    get_vertices(q).currentk += -1
    return
end

@inline function get_indices(q::QMC)
    return get_vertices(q).indices
end

@inline function get_values(q::QMC)
    return get_vertices(q).values
end

@inline function get_expgamma0(q::QMC, i)
    return @inbounds q.p.expgamma0[i]
end

@inline function get_L_temp(q::QMC, k)
    return view(q.v.L_temp, 1:k)
end

@inline function get_R_temp(q::QMC, k)
    return view(q.v.R_temp, 1:k)
end

@inline function get_Mklm_temp(q::QMC, k)
    return view(q.v.Mklm_temp, 1:k, 1:k, :)
end

@inline function get_ev_temp(q::QMC)
    return q.v.ev_temp
end

@inline function get_g0l_temp(q::QMC, k)
    return view(q.v.g0l_temp, 1:k, :)
end

@inline function get_ntemp_temp(q::QMC, k)
    return view(q.v.ntemp_temp, 1:k, 1:k)
end

@inline function get_Mkl_temp(q::QMC)
    return q.v.Mkl_temp
end


function MP_test(ratio)
    r = rand()
    pass = ifelse(min(1.0, abs(ratio)) > r, true, false)
    return pass
end

@inline function updateDc!(q, value, i, ispin)
    index = get_indices(q)[i]
    @inbounds q.v.Dc[index, ispin] = value
    return
end

@inline function updateDr!(q, value, i, ispin)
    index = get_indices(q)[i]
    q.v.Dr[index, ispin] = value
    return
end

@inline function updateN_direct!(q, values, ispin)
    n = size(values)
    @inbounds for j = 1:n[2]
        for i = 1:n[1]
            q.v.N[i, j, ispin] = values[i, j]
        end
    end
    #q.v.N[1:n[1],1:n[2],ispin] = values[1:n[1],1:n[2]]
    return
end

@inline function updateN_direct!(q, value::Float64, i, j, ispin)
    @inbounds q.v.N[i, j, ispin] = value
    return
end

function dot_Dr_N_Dc(q::QMC, k, ispin)
    Dr = get_Dr(q, k, ispin)
    Dc = get_Dc(q, k, ispin)
    N = view(get_N(q), 1:k, 1:k, ispin)
    temp = get_R_temp(q, k)
    mul!(temp, N, Dc)
    return dot(Dr, temp)

    #return dot(view(q.v.Dr,1:k,ispin),view(q.v.N,1:k,1:k,ispin)*view(q.v.Dc,1:k,ispin))
end

@inline function calc_Gf0_tau(q::QMC, tau, l)
    return calc_Gf0_tau(q.p.gtau0, tau, l)
end

function calc_Gf0_tau(q::QMC, tau_i, tau_j, l)
    dtau = tau_i - tau_j
    if dtau < 0
        dtau += get_beta(q)
        Gtau0 = -calc_Gf0_tau(q, dtau, l)
    elseif dtau == 0.0
        Gtau0 = calc_Gf0_tau(q, 1e-8, l)
    else
        Gtau0 = calc_Gf0_tau(q, dtau, l)
    end
    return Gtau0
end

function generate_remove_vertex(q::QMC)
    index = rand(1:get_currentk(q))
    vertex = get_values(q)[get_indices(q)[index]]
    return vertex, index
end

function insert_vertex!(q::QMC, vertex, position)
    insert_vertex!(q.v.vertices, vertex, position)
    return
end

function remove_vertex!(q::QMC, vertex, position)
    remove_vertex!(q.v.vertices, vertex, position)
    return
end

function remove_vertex!(v::Vertices, vertex, position)
    is = v.indices[position]
    deleteat!(v.values, is)
    deleteat!(v.indices, position)
    for i in 1:v.currentk-1
        v.indices[i] += ifelse(v.indices[i] >= is, -1, 0)
    end

    return
end

function remove_vertex!(v::Vertices, tau, spin, x, position)
    remove_vertex!(v, Vertex(tau, spin, x), position)
    return
end


function calc_DσΓσ(q::QMC, ispin)
    k = get_currentk(q)
    Dσ = zeros(k, k)
    Γσ = zeros(k, k)
    Ninv = calc_Ninv(q, ispin)
    for i = 1:k
        Dσ[i, i] = Ninv[i, i]
        for j = 1:k
            if i != j
                Γσ[j, i] = -Ninv[j, i]
            end
        end
    end
    return Dσ, Γσ
end

function calc_logdet!(q::QMC)
    ispin = 1
    detw_up = calc_logdet(q, ispin)
    ispin = 2
    detw_down = calc_logdet(q, ispin)
    logdetW = log(detw_up * detw_down)
    return logdetW
end

function calc_logdet(q::QMC, ispin)
    Ninv = calc_Ninv(q, ispin)
    return det(Ninv)
end


function calc_Ninv(q::QMC, ispin)
    k = get_currentk(q)
    ssign = (-1)^(ispin - 1)
    exp_v = zeros(Float64, k)
    indices = get_indices(q)
    for j = 1:k
        vertex_j = get_timeordered_vertex(q, j)
        s_j = vertex_j.spin
        exp_v[indices[j]] = ifelse(ssign * s_j == 1, get_expgamma0(q, 1), get_expgamma0(q, 2))
    end
    Gmm = zeros(k, k)
    Ninv = zeros(k, k)

    for j = 1:k
        jj = get_indices(q)[j]
        tau_j = get_timeordered_vertex(q, j).tau
        for i = 1:k
            ii = get_indices(q)[i]
            tau_i = get_timeordered_vertex(q, i).tau
            #dtau = tau_i - tau_j
            Gmm[ii, jj] = calc_Gf0_tau(q, tau_i, tau_j, ispin)
        end
    end

    for j = 1:k
        #tau_j = get_timeordered_vertex(q, j).tau
        #tau_j = get_timeordered_vertex(q, j).tau
        for i = 1:k
            #tau_i = get_timeordered_vertex(q, i).tau
            if i == j
                Ninv[i, j] += exp_v[i]
            end
            Ninv[i, j] += -Gmm[i, j] * (exp_v[j] - 1)
        end
    end
    return Ninv
end
