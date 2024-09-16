function I0(P,h1,h2,h3,h4)
    #calculate 1D PBF

    hh = h1*h1 + h2*h2 + h3*h3 + h4*h4
    R = sqrt(P*P + hh)
    hh1 = hh - h1*h1
    zero2 = zerotol*zerotol

    # case = cases(h1,h2,h3,h4)
    # print("case: $case")
    if testMode == 1
        println("hs: ", h1, " ", h2, " ", h3, " ", h4, " ")
    end
    
    if hh < zero2*P*P
        Phi1 = log(P) / P
        I = Phi1 / 6
    else
        Phi1 = 0.5 * log((P + R) * (P + R) / hh) / P
        if hh1 < zero2 * hh
            I = Phi1 / 6
        else
            if h3*h3 + h4*h4 < zero2 * hh
                if h1*P == 0
                    Phi2 = 1 / (hh + h2*R)
                else
                    Phi2 = atan(h1*P / (hh + h2*R)) / (h1*P)
                end
                I = (Phi1 - h2*Phi2) / 6
            else
                R1 = sqrt(P*P + h1*h1)
                if h2*h2 + h4*h4 < zero2 * hh
                    if h1*P == 0
                        Phi2 = 1 / (hh + h3*R)
                    else
                        Phi2 = atan(h1*P / (hh + h3*R)) / (h1*P)
                    end
                    Phi3 = 0.5 * hh1 / R1 * log((R1 + R) * (R1 + R) / hh1)
                    I = ((h1*h1 - h3*h3) * Phi1 - 2*h1*h1*h3*Phi2 + Phi3) / (6 * h1*h1)
                elseif h2*h2 + h3*h3 < zero2 * hh
                    if h1*P == 0
                        Phi2 = 1 / (hh + h4*R)
                    else
                        Phi2 = atan(h1*P / (hh + h4*R)) / (h1*P)
                    end
                    Phi3 = 0.5 / R1 * log((R1 + R) * (R1 + R) / hh1)
                    I = ((h1*h1 - 3*h4*h4) * Phi1 - h4 * (3*h1*h1 - h4*h4) * Phi2 + 3*h4*h4*Phi3 - h4*h4 / (R + h4)) / (6 * h1*h1)
                elseif h1*h1 + h4*h4 < zero2 * hh
                    h = sqrt(hh)
                    R2 = sqrt(P*P + h2*h2)
                    Phi4 = h3*h3 / (h2 * P*P) * ((R2 / h2 * log((R2 + R) / h3) - log((h2 + h) / h3)))
                    I = (hh / (h2*h2) * Phi1 - 1 / (R + h) - Phi4) / 6
                elseif h1*h1 + h3*h3 < zero2 * hh
                    h = sqrt(hh)
                    R2 = sqrt(P*P + h2*h2)
                    Phi2 = atan(h2*P / (hh + h4*R)) / P
                    Phi4 = h4*h4 / (h2 * P*P) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))
                    I = ((1 + 3*h4*h4 / (h2*h2)) * Phi1 - 2 * (h4 / h2)^3 * Phi2 - 3 * Phi4 + (2*h4*h4 - h2*h2) / (h2*h2 * (R + h))) / 6
                else
                    I = 0
                    if testMode == 1
                    println("I0 executed a prohibited option")
                    println("h1 = $(h1)  h2 = $(h2)  h3 = $(h3)  h4 = $(h4)")
                    end
                end
            end
        end
    end

    if testMode == 1
        println("I0: ", I)
    end

    return I

end