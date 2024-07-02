function Identify_cases(h1, h2, h3, h4)
    if h1 > zerotol && h2 < zerotol && h3 < zerotol && h4 < zerotol
        global case = 1
    elseif h1 < zerotol && h2 > zerotol && h3 < zerotol && h4 < zerotol
        global case = 2
    elseif h1 > zerotol && h2 > zerotol && h3 < zerotol && h4 < zerotol
        global case = 3
    elseif h1 > zerotol && h2 < zerotol && h3 > zerotol && h4 < zerotol
        global case = 4
    elseif h1 < zerotol && h2 > zerotol && h3 > zerotol && h4 < zerotol
        global case = 5
    elseif h1 > zerotol && h2 < zerotol && h3 < zerotol && h4 > zerotol
        global case = 6
    elseif h1 < zerotol && h2 > zerotol && h3 < zerotol && h4 > zerotol
        global case = 7
    else
        global case = 8 #case 8 sould give the error below
    end
    println("The case is $case")
    if case >= 8
        error("Error: h does not satisfy any case")
    end
    return case
end