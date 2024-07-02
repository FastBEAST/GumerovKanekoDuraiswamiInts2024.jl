function I1(e1, a1, h2, h3, h4) :: Float64
    #evaluation of 1D integral
    #the a passed down here should have 1 vector in the format of Vector{Vector{Float64}}

    T = eltype(e1)
    #find parameters h1 and s01
    h1, s01 = GSorthogonalization_expan(e1,a1)
    s01 = s01[1]

    #making the smallest two h zero
    h_all = [h1,h2,h3,h4]
    i = sortperm(h_all)
    h_all[i[1]] = T(0); h_all[i[2]] = T(0)
    h1 = h_all[1]; h2 = h_all[2]; h3 = h_all[3]; h4 = h_all[4]

    na1 = norm(a1)

    I = T(0)
    P1 = abs(T(1)+s01)*na1
    P2 = abs(s01)*na1
    if abs(T(1) + s01) > zerotol
        I += (T(1) + s01) * I0(P1, h1, h2, h3, h4)
    end
    if abs(s01) > zerotol
        I -= s01 * I0(P2, h1, h2, h3, h4)
    end
    
    if testMode == 1
        println("I1: $I")
    end

    return I

end