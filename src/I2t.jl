function I2t(e, a, h3, h4)
    #Find 2D integral from 2D domian boundaries(1D integrals)
    #There are triangle surfaces and square surfaces. This is for the triangle surfaces.
    #The a passed down here should have 2 vectors.

    T = eltype(e)
    #find parameters h2 and s0
    h2, s0 = GSorthogonalization_expan(e, a)

    #Define parmeters of the 1D integrals
    a12 = a[2]
    a14 = a[1]
    a15 = a[1] - a[2]
    e12 = e14 = a[1] * s0[1] + a[2] * s0[2]
    e15 = e12 + a[2]


    #calculate boundary integrals
    if abs(s0[1]) < zerotol
        I11t = T(0)
    else
        I11t = -s0[1] * I1(e12, SVector{1,Vector{Float64}}([a12]), h2, h3, h4)
    end
    if abs(s0[2]) < zerotol
        I12t = T(0)
    else
        I12t = -s0[2] * I1(e14, SVector{1,Vector{Float64}}([a14]), h2, h3, h4)
    end
    if abs(T(1) + s0[1] + s0[2]) < zerotol
        I13t = T(0)
    else
        I13t = (T(1) + s0[1] + s0[2]) * I1(e15, SVector{1,Vector{Float64}}([a15]), h2, h3, h4)
    end

    #calculate 2D integral for square
    I2t_sum = I11t + I12t + I13t

    return I2t_sum
end
