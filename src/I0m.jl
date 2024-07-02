function I0m(P,h1,h2,h3,h4)
    # calculate 1D PBF

    zero2 = zerotol^2
    hh = h1^2 + h2^2 + h4^2
    R = sqrt(P^2 + hh)
    Phi1 = 0.5 * log((P + R) * (P + R) / hh) / P
    R1 = sqrt(P*P + h1*h1)

    if testMode == 1
        println("hsm: ", h1, " ", h2, " ", h3, " ", h4, " ")
    end
    
    if h2*h2 < zero2 * hh
        #case 6
        if P == 0.0 #< zerotol??
            Phi2 = h1 / (hh + h4 * R)
        else
            Phi2 = atan(h1 * P / (hh + h4 * R)) / P
        end
        Phi3 = 1 / R1 * log((R1 + R) / h4)
        I = (Phi1 + (h1*h1 - h4*h4) / (2 * h1 * h4) * Phi2 - Phi3 + 0.5 / (R + h4)) / (h1*h1)
    elseif h1*h1 < zero2 * hh
        #case 7
        h = sqrt(hh)
        R2 = sqrt(P*P + h2*h2)
        Phi2 = atan(h2 * P / (hh + h4 * R)) / P
        Phi4 = h4*h4 / (h2 * P*P) * ((R2 / h2 * log((R2 + R) / h4) - log((h2 + h) / h4)))
        I = (-Phi1 + h4 / h2 * Phi2 + h2*h2 / (h4*h4) * Phi4 - (R2*R2 / (R + h4) - h2*h2 / (h + h4)) / (P*P)) / (h2*h2)
    else
        I = 0
        if testMode == 1
        println("I0 executed a prohibited option")
        println("h1 = $h1  h2 = $h2  h3 = $h3  h4 = $h4")
        end
        # error("Execution stopped")
    end
    
    if (isnan(I) || isinf(I)) && testMode == 1
        println("NaN or inf")
    end
    
    return I
end