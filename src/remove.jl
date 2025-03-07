function qmc_remove!(q::QMC)
    if get_currentk(q) == 0
        pass = false
        rsign = 1.0
        return pass, rsign
    end
    vertex, position = generate_remove_vertex(q)
    det_ratio = calc_remove_ratio!(q, vertex, position)
    ratio0 = det_ratio[1] * det_ratio[2]
    ratio = (get_currentk(q) / get_K(q)) * ratio0
    rsign = sign(ratio)
    #println("ratio_remove: $ratio")

    pass = MP_test(ratio)
    if pass
        index = get_indices(q)[position]
        remove_vertex!(q, vertex, position)
        currentk_remove1!(q)
        update_remove_N!(q, index, det_ratio)
        q.v.P[get_currentk(q)+1] += 1
        q.v.current_logW += log(ratio0)
    end

    return pass, rsign
end


function update_remove_N!(q::QMC, position, det_ratio)
    update_remove_N!(q, position, 1, det_ratio[1])
    update_remove_N!(q, position, 2, det_ratio[2])
    return
end

function update_remove_N!(q::QMC, index, rspin, det_ratio)
    #        println("position: $position")
    as = index
    k = get_currentk(q)
    if k == 0
        return
    end


    lambda = det_ratio
    ntemp = get_ntemp_temp(q, k) #zeros(Float64,k,k)
    N = get_N(q)
    @inbounds for j in 1:k
        jj = j + ifelse(j - as >= 0, 1, 0)
        for i in 1:k
            ii = i + ifelse(i - as >= 0, 1, 0)
            ntemp[i, j] = N[ii, jj, rspin] - N[ii, as, rspin] * N[as, jj, rspin] / lambda
        end
    end

    updateN_direct!(q, ntemp, rspin)
    return
end

function calc_remove_ratio!(q, vertex, position)
    index = get_indices(q)[position]
    det_ratio_up = get_N(q)[index, index, 1]
    det_ratio_down = get_N(q)[index, index, 2]
    return det_ratio_up, det_ratio_down
end
