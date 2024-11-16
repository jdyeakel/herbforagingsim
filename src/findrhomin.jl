function findrhomin(rhoexpvec,probsurv)
    
    rhovec = 10 .^ rhoexpvec;

    first_nonzero_pos = findfirst(x -> x > 0.0, probsurv)
    low_surv_pos = findall(x -> (x > 0.0 && x < 0.5), probsurv)


    # Retrieve the corresponding value from `rhovec`
    rhomin_values = first_nonzero_pos !== nothing ? rhovec[low_surv_pos] : maximum(rhovec) #missing
    mean_rhomin = mean(rhomin_values)

    return mean_rhomin

end