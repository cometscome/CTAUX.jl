function make_hybridization(mesh_omega, V, norb)
    mfreq = length(mesh_omega)
    delta = zeros(ComplexF64, mfreq, norb, norb)
    for i = 1:mfreq
        z = mesh_omega[i] * im
        for j = 1:norb
            delta[i, j, j] = 1
        end
        delta[i, :, :] = delta[i, :, :] * (z - im * sqrt(1 - z^2)) / 2
    end
    delta = delta * V^2
    return delta
end

function make_Gf0_omega(mesh_omega, norb, Ef, delta, U)
    mfreq = length(mesh_omega)
    Gf0 = zeros(ComplexF64, mfreq, norb, norb) #Non-perturbative Green's function in the Matsubara space
    for i in 1:mfreq
        for j in 1:norb
            Gf0[i, j, j] = 1 / (im * mesh_omega[i] - Ef - delta[i, j, j] - U / 2)
        end
    end

    return Gf0
end

function make_Gf0_tau(Gf0_omega, mesh_omega, ntime, beta)
    return fft_ω2τ(Gf0_omega, mesh_omega, ntime, beta)
end


function set_greenfunctions(Gf0_tau, mesh_tau)
    nl = size(Gf0_tau)[2]
    println(nl)
    gtau0 = Greenfunction[]
    #        gtau0 = Array{Greenfunction}{undef,nl}
    for l = 1:nl
        spl = Dierckx.Spline1D(mesh_tau, Gf0_tau[:, l, l])
        push!(gtau0, Greenfunction(spl))
    end
    return gtau0

end

function fft_ω2τ(f_omega, mesh_omega, ntime, beta)
    mfreq = size(f_omega)[1]
    nl = size(f_omega)[2]
    #        println("nl = $nl")
    f_tau = zeros(Float64, ntime, nl, nl)
    for i = 1:nl
        for j = 1:nl
            f_tau[:, i, j] = fft_backward(ntime, beta, f_omega[:, i, j], mesh_omega)
        end
    end

    for i = 1:nl
        for itau = 1:ntime
            if f_tau[itau, i, i] > 0
                f_tau[itau, i, i] = 1e-6 #to avoid positive values
            end
        end
    end
    return f_tau
end

function fft_backward(ntime, beta, v_omega, mesh_omega)
    mfreq = length(mesh_omega)
    tail = calc_tails(v_omega, mesh_omega)
    gk = zeros(ComplexF64, 2mfreq)

    for j = 1:mfreq
        gk[2j] = v_omega[j] - tail / (im * mesh_omega[j])
    end

    fft!(gk)


    g_tau = zeros(Float64, ntime)
    g_tau[1:ntime] = real(gk[1:ntime]) * (2 / beta) .- tail / 2

    a = real(v_omega[mfreq]) * mesh_omega[mfreq] / π
    g_tau[1] += a
    g_tau[ntime] += -a

    return g_tau

end

function calc_tails(vω, mesh_omega) #to calculate a tail
    mfreq = length(mesh_omega)
    ntail = 128

    Sn = 0.0
    Sx = 0.0
    Sy = 0.0

    Sxx = 0.0
    Sxy = 0.0

    for j in mfreq-ntail:mfreq
        ωn = mesh_omega[j]
        Sn += 1
        Sx += 1 / ωn^2
        Sy += imag(vω[j]) * ωn
        Sxx += 1 / ωn^4
        Sxy += imag(vω[j]) * ωn / ωn^2
    end
    rtail = (Sx * Sxy - Sxx * Sy) / (Sn * Sxx - Sx * Sx)

    return rtail
end


@inline function calc_Gf0_tau(gtau0::Array{Greenfunction,1}, tau, l)
    return gtau0[l].spl(tau)
end

function calc_green(S, τmesh, β, ntime, norbs, Gτ0spl)
    #        global β,ntime,norbs
    δτ = β / (ntime - 1)
    Gτ = zeros(Float64, ntime, norbs)
    gsum = zeros(Float64, norbs)
    g0 = zeros(Float64, norbs)


    for i in 1:ntime
        τ = τmesh[i]
        if i == ntime
            τ += -1e8
        elseif i == 1
            τ += 1e-8
        end

        for sigma in 1:norbs
            Gτ[i, sigma] = Gτ0spl[sigma](τ)
        end
        gsum[:] .= 0.0
        for j in 1:i
            τt = τmesh[j]
            dτ = τ - τt
            if i == j
                dτ = 1e-8
            end

            for sigma in 1:2
                g0[sigma] = Gτ0spl[sigma](dτ)
                gsum[sigma] += g0[sigma] * S[j, sigma] * δτ * ifelse(j == 1 || j == i, 0.5, 1.0)

                #=
                if j == 1 || j ==i
                    gsum[sigma] += g0[sigma]*S[j,sigma]*δτ/2
                else
                    gsum[sigma] += g0[sigma]*S[j,sigma]*δτ
                end
                =#


            end

        end
        for j in i:ntime
            τt = τmesh[j]
            dτ = τ - τt
            if i == j
                dτ = -1e-8
            end
            for sigma in 1:2
                g0[sigma] = -Gτ0spl[sigma](dτ + β)
                gsum[sigma] += g0[sigma] * S[j, sigma] * δτ * ifelse(j == ntime || j == i, 0.5, 1.0)
                #=
                if j == ntime || j ==i
                    gsum[sigma] += g0[sigma]*S[j,sigma]*δτ/2
                else
                    gsum[sigma] += g0[sigma]*S[j,sigma]*δτ
                end
                =#


            end


        end
        for sigma in 1:norbs
            Gτ[i, sigma] += gsum[sigma]
        end
    end


    return Gτ
end

