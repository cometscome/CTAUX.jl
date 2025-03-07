function calc_linear_mesh(xmin, xmax, n)
    x = zeros(Float64, n)
    for i in 1:n
        x[i] = (xmax - xmin) * (i - 1) / (n - 1) + xmin
    end
    return x
end