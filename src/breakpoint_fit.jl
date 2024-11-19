function breakpoint_fit(massvec,rhomin,break_mass,index_break,zetaindex)
    
    massvec1 = massvec[1:index_break]
    rhomin1 = rhomin[1:index_break,zetaindex]
    massvec2 = massvec[index_break+1:end]
    rhomin2 = rhomin[index_break+1:end,zetaindex]

    log_massvec1 = log.(massvec1)
    log_rhomin1 = log.(rhomin1)
    # Fit a line
    X1 = hcat(ones(length(log_massvec1)), log_massvec1)
    coeffs1 = X1 \ log_rhomin1  # Linear regression
    log_a1, b1 = coeffs1
    a1 = exp(log_a1)  # Convert log(a1) to a1

    log_massvec2 = log.(massvec2)
    log_rhomin2 = log.(rhomin2)
    # Fit a line
    X2 = hcat(ones(length(log_massvec2)), log_massvec2)
    coeffs2 = X2 \ log_rhomin2  # Linear regression
    log_a2, b2 = coeffs2
    a2 = exp(log_a2)  # Convert log(a2) to a2

    mass_range1 = range(minimum(massvec), stop=break_mass, length=100)
    rhomin_fit1 = a1 .* mass_range1 .^ b1;

    mass_range2 = range(break_mass, stop= maximum(massvec))
    rhomin_fit2 = a2 .* mass_range2 .^ b2;

    return a1, b1, a2, b2, mass_range1, rhomin_fit1, mass_range2, rhomin_fit2

end