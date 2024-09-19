function GalerkinLaplaceTriGS(x1, x2, x3, y1, y2, y3; L_output=1, Md_output=0)
    #L means single layer, M double layer, Ld normal derivative of single layer (transposed double layer), Md normal derivative for double layer (hypersingular)
    #set L_output to 1 to calculate L, M, Ld; set Md_output to 1 to calculate Md; set them to 0 to output zero. This saves time if not all operators are needed.
    #Md has not been tested with BEAST.jl

    T = promote_type(eltype(x1), eltype(x2), eltype(x3), eltype(y1), eltype(y2), eltype(y3))

    #Triangle transformation
    # Define the edge vectors
    lx1 = x2 - x1
    lx2 = x3 - x2
    lx3 = x1 - x3
    ly1 = y2 - y1
    ly2 = y3 - y2
    ly3 = y1 - y3

    # Calculate the normal vectors and areas
    nx = cross(lx1, -lx3)
    Ax = norm(nx)
    nx = nx / Ax
    Ax = T(0.5) * Ax
    ny = cross(ly1, -ly3)
    Ay = norm(ny)
    ny = ny / Ay
    Ay = T(0.5) * Ay

    # Calculate the lengths of the triangle's edges
    ellx1 = norm(lx1)
    ellx2 = norm(lx2)
    ellx3 = norm(lx3)
    elly1 = norm(ly1)
    elly2 = norm(ly2)
    elly3 = norm(ly3)

    # Calculate the unit normal vectors
    ncx1 = cross(lx1, nx)
    ncx2 = cross(lx2, nx)
    ncx3 = cross(lx3, nx)
    ncx1 = ncx1 / norm(ncx1)
    ncx2 = ncx2 / norm(ncx2)
    ncx3 = ncx3 / norm(ncx3)
    ncy1 = cross(ly1, ny)
    ncy2 = cross(ly2, ny)
    ncy3 = cross(ly3, ny)
    ncy1 = ncy1 / norm(ncy1)
    ncy2 = ncy2 / norm(ncy2)
    ncy3 = ncy3 / norm(ncy3)

    # Calculate the dot product
    nxy = dot(nx, ny)
    nxy1 = min(norm(nx + ny), norm(nx - ny))

    a1 = lx1
    a2 = -lx3
    a3 = -ly1
    a4 = ly3
    e4 = x1 - y1
    h4, s0 = GSorthogonalization_expan(e4, SVector(a1, a2, a3, a4))

    #Define parmeters of the 3D integrals
    a11 = a12 = a13 = a34 = a1
    a21 = a22 = a23 = a36 = a2
    a31 = a14 = a15 = a16 = a3
    a33 = a24 = a25 = a26 = a4
    a32 = a4 - a3
    a35 = a2 - a1
    e31 = e33 = e34 = e36 = a1 * s0[1] + a2 * s0[2] + a3 * s0[3] + a4 * s0[4]
    e32 = e31 + a3
    e35 = e31 + a1

    if L_output == 1
        #calculate single layer
        I31 = I3(e31, SVector(a11, a21, a31), h4)
        I32 = I3(e32, SVector(a12, a22, a32), h4)
        I33 = I3(e33, SVector(a13, a23, a33), h4)
        I34 = I3(e34, SVector(a14, a24, a34), h4)
        I35 = I3(e35, SVector(a15, a25, a35), h4)
        I36 = I3(e36, SVector(a16, a26, a36), h4)
        I4 = -s0[4] * I31 + (T(1) + s0[3] + s0[4]) * I32 - s0[3] * I33 - s0[2] * I34 + (T(1) + s0[1] + s0[2]) * I35 - s0[1] * I36
        L = T(4) * Ax * Ay * I4 #jacobian is 2A for transformation so here 2Ax * 2Ay = 4 * Ax * Ay

        #caclulate double layer
        if nxy1 > 1e-9
            #nondegenerate case
            Fx = 6 * Ay * (ellx1 * ncx1 * I34 + ellx3 * ncx3 * I36 + ellx2 * ncx2 * I35) # I'= 3I from the paper formula (4.8), then use (4.5)
            Fy = 6 * Ax * (elly1 * ncy1 * I31 + elly3 * ncy3 * I33 + elly2 * ncy2 * I32)

            M = (dot(nx, Fy) - nxy * dot(ny, Fx)) / (nxy^2 - 1) #may not work but worked fine till now
        else
            #degenerate case
            delta = -dot(e4, nx)#singed distance between planes

            if h4 <= zerotol
                M = T(0.0)
            else
                I31m = I3m(e31, [a11, a21, a31], h4)
                I32m = I3m(e32, SVector(a12, a22, a32), h4)
                I33m = I3m(e33, [a13, a23, a33], h4)
                I34m = I3m(e34, SVector(a14, a24, a34), h4)
                I35m = I3m(e35, SVector(a15, a25, a35), h4)
                I36m = I3m(e36, SVector(a16, a26, a36), h4)

                I4m =
                    -s0[4] * I31m + (1 + s0[3] + s0[4]) * I32m - s0[3] * I33m - s0[2] * I34m + (1 + s0[1] + s0[2]) * I35m -
                    s0[1] * I36m
                M = T(4) * Ax * Ay * delta * I4m
            end


            I34d = I3(e4, SVector(a14, a24, a34), T(0))
            I35d = I3(e4 + a1, SVector(a15, a25, a35), T(0))
            I36d = I3(e4, SVector(a16, a26, a36), T(0))
            Fx = 6 * Ay * (ellx1 * ncx1 * I34d + ellx3 * ncx3 * I36d + ellx2 * ncx2 * I35d)
        end

        #double layer transposed
        Ld = (-Fx - M * nx) ⋅ (ny)
    else
        L = 0.0
        M = 0.0
        Ld = 0.0
    end

    #finding normal derivative for double layer
    if Md_output == 1
        H11 = I2s(x1 - y1, SVector(lx1, -ly1), T(0), T(0))
        H12 = I2s(x1 - y2, SVector(lx1, -ly2), T(0), T(0))
        H13 = I2s(x1 - y3, SVector(lx1, -ly3), T(0), T(0))
        H21 = I2s(x2 - y1, SVector(lx2, -ly1), T(0), T(0))
        H22 = I2s(x2 - y2, SVector(lx2, -ly2), T(0), T(0))
        H23 = I2s(x2 - y3, SVector(lx2, -ly3), T(0), T(0))
        H31 = I2s(x3 - y1, SVector(lx3, -ly1), T(0), T(0))
        H32 = I2s(x3 - y2, SVector(lx3, -ly2), T(0), T(0))
        H33 = I2s(x3 - y3, SVector(lx3, -ly3), T(0), T(0))

        Md =
            -6 * (
                dot(lx1, ly1) * H11 +
                dot(lx1, ly2) * H12 +
                dot(lx1, ly3) * H13 +
                dot(lx2, ly1) * H21 +
                dot(lx2, ly2) * H22 +
                dot(lx2, ly3) * H23 +
                dot(lx3, ly1) * H31 +
                dot(lx3, ly2) * H32 +
                dot(lx3, ly3) * H33
            )
    else
        Md = 0.0
    end

    return L, M, Ld, Md

end
