function I3m(e, a, h4)
    #Find 3D integral from 3D domian boundaries(2D integrals)
    #the 3D shapes are prisms; The a passed down here should have 3 vectors

    T = eltype(e)
    #find h3 and s0
    h3, s0 = GSorthogonalization_expan(e, a)

    #Define parmeters of the 2D integrals
    a11 = a[1] - a[2]
    a21 = a22 = a23 = a[3]
    a12 = a24 = a25 = a[2]
    a13 = a14 = a15 = a[1]
    e22 = e23 = e25 = a[1] * s0[1] + a[2] * s0[2] + a[3] * s0[3]
    e21 = e22 + a[2]
    e24 = e22 + a[3]

    #calculate boundary integrals
    if abs(1 + s0[1] + s0[2]) < zerotol
        I21 = T(0)
    else
        I21 = (T(1) + s0[1] + s0[2]) * I2sm(e21, SVector(a11, a21), h3, h4)
    end
    if abs(s0[1]) < zerotol
        I22 = T(0)
    else
        I22 = -s0[1] * I2sm(e22, SVector(a12, a22), h3, h4)
    end
    if abs(s0[2]) < zerotol
        I23 = T(0)
    else
        I23 = -s0[2] * I2sm(e23, SVector(a13, a23), h3, h4)
    end
    if abs(T(1) + s0[3]) < zerotol
        I24 = T(0)
    else
        I24 = (T(1) + s0[3]) * I2tm(e24, SVector(a14, a24), h3, h4)
    end
    if abs(s0[3]) < zerotol
        I25 = T(0)
    else
        I25 = -s0[3] * I2tm(e25, SVector(a15, a25), h3, h4)
    end

    #calculate 3D integral
    I3_sum = I21 + I22 + I23 + I24 + I25

    return I3_sum

end
