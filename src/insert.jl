function qmc_insert!(q::QMC)
    vertex, position = generate_insert_vertex(q)
    det_ratio = calc_insert_ratio!(q, vertex)
    #println(det_ratio)
    ratio = det_ratio[1] * det_ratio[2] * get_K(q) / (get_currentk(q) + 1)

    #println("det_ratioup and down","$det_ratio")
    rsign = sign(ratio)
    pass = MP_test(ratio)
    #println("ratio: $ratio")

    if pass
        insert_vertex!(q, vertex, position)
        currentk_add1!(q)
        update_insert_N!(q, position, det_ratio)

        q.v.P[get_currentk(q)+1] += 1

    end

    return pass, rsign
end



function generate_insert_vertex(q::QMC)
    tau = rand() * get_beta(q)
    spin = (-1)^(rand(1:2) - 1)
    #        x = rand(1:get_L(q))
    vertex = Vertex(tau, spin)

    index = 1
    currentk = get_currentk(q)
    indices = get_indices(q)
    values = get_values(q)

    if currentk > 0
        if tau < values[indices[1]].tau
            index = 1
        elseif tau > values[indices[currentk]].tau #τ > τcon[indexcon[currentk]]
            index = currentk + 1
        else
            i = 1
            for j in 1:currentk
                i += ifelse(values[indices[j]].tau < tau, 1, 0)
            end
            index = i
        end
    end

    return vertex, index
end

function calc_insert_ratio!(q::QMC, vertex)
    ispin = 1
    det_ratio_up = calc_insert_ratio_spin!(q, vertex, ispin)
    ispin = 2
    det_ratio_down = calc_insert_ratio_spin!(q, vertex, ispin)
    return det_ratio_up, det_ratio_down
end

function calc_insert_ratio_spin!(q::QMC, vertex, ispin)
    k = get_currentk(q)
    ssign = (-1)^(ispin - 1)
    #println(get_expgamma0(q,1))
    ev = ifelse(ssign * vertex.spin == 1, get_expgamma0(q, 1), get_expgamma0(q, 2))

    Gtau0 = calc_Gf0_tau(q, 1e-8, ispin)
    #println("Gtau0 $Gtau0") 
    Dd = ev - Gtau0 * (ev - 1)
    #println("Dd ", Dd)

    if k == 0
        det_ratio = Dd
        return det_ratio
    end

    tau_j = vertex.tau
    sj = vertex.spin
    ev = ifelse(ssign * sj == 1, get_expgamma0(q, 1), get_expgamma0(q, 2))
    for i in 1:k
        tau_i = get_timeordered_vertex(q, i).tau
        Gtau0 = calc_Gf0_tau(q, tau_i, tau_j, ispin)
        value = -Gtau0 * (ev - 1)
        updateDc!(q, value, i, ispin)
    end

    tau_i = vertex.tau
    for j = 1:k
        vertex_j = get_timeordered_vertex(q, j)
        tau_j = vertex_j.tau
        s_j = vertex_j.spin
        Gtau0 = calc_Gf0_tau(q, tau_i, tau_j, ispin)
        ev = ifelse(ssign * s_j == 1, get_expgamma0(q, 1), get_expgamma0(q, 2))
        value = -Gtau0 * (ev - 1)
        updateDr!(q, value, j, ispin)
    end
    lambda = Dd - dot_Dr_N_Dc(q, k, ispin)

    det_ratio = lambda[1]

    return det_ratio
end


function update_insert_N!(q::QMC, position, det_ratio)
    update_insert_N!(q, position, 1, det_ratio[1])
    update_insert_N!(q, position, 2, det_ratio[2])
    return
end

function update_insert_N!(q::QMC, position, rspin, det_ratio)
    k = get_currentk(q)
    if k == 1
        updateN_direct!(q, 1 / det_ratio, 1, 1, rspin)
        return
    end
    lambda = det_ratio

    #R = zeros(Float64,k-1)
    R = get_R_temp(q, k - 1)
    N = view(get_N(q), 1:k-1, 1:k-1, rspin)
    Dr = view(get_Dr(q), 1:k-1, rspin)
    mul!(R, N', Dr)
    R .*= -1 / lambda
    #R = -N'*get_Dr(q)[1:k-1,rspin]/lambda
    #R = - nmat[1:k-1,1:k-1,rspin]'*Dr[1:k-1,rspin]/λ

    #L = zeros(Float64,k-1)
    L = get_L_temp(q, k - 1)
    Dc = view(get_Dc(q), 1:k-1, rspin)
    mul!(L, N, Dc)
    L .*= -1 / lambda
    #L = - N*get_Dc(q)[1:k-1,rspin]/lambda
    #        nmat[1:k-1,1:k-1,rspin]*Dc[1:k-1,rspin]/λ

    ntemp = get_ntemp_temp(q, k) #zeros(Float64,k,k)
    #nmatt = zeros(Float64,k-1,k-1)
    #nmatt=N

    for j in 1:k-1
        for i in 1:k-1
            N[i, j] += lambda * L[i] * R[j]
        end
    end
    as = get_indices(q)[position] #indexcon[vertex[3]]        
    ntemp[as, as] = 1 / lambda
    #println("as ",as,"\t",k)
    ntemp[1:k-1, 1:k-1] .= view(N, 1:k-1, 1:k-1)
    ntemp[k, 1:k-1] .= view(R, 1:k-1)
    ntemp[1:k-1, k] .= view(L, 1:k-1)

    updateN_direct!(q, ntemp, rspin)

    #        nmat[1:k,1:k,rspin] = ntemp[1:k,1:k]

    return
end


function insert_vertex!(v::Vertices, tau, spin, position)
    insert_vertex!(v, Vertex(tau, spin), position)

    return
end

function insert_vertex!(v::Vertices, vertex, position)
    if v.currentk == 0
        index = 1
        insert!(v.values, 1, vertex)
        insert!(v.indices, 1, index)
        return
    end
    index = v.currentk + 1
    insert!(v.indices, position, index)       #insert index at position 
    insert!(v.values, index, vertex)          #insert vertex at index
    return
end
