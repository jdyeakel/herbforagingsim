function findrhomin(rhoexpvec,probsurv,surv_threshold)
    
    rhovec = 10 .^ rhoexpvec;

    first_nonzero_pos = findfirst(x -> x > surv_threshold, probsurv)
    rhomin = first_nonzero_pos !== nothing ? rhovec[first_nonzero_pos] : maximum(rhovec) #missing
    return rhomin

    # #NOTE: Need to think about this!
    # #Issues bc a lot of probsurv vectors go from 0 -> 0.8, etc.
    # low_surv_pos = findall(x -> (x > 0.0 && x < 0.5), probsurv)
    # # Retrieve the corresponding value from `rhovec`
    # rhomin_values = first_nonzero_pos !== nothing ? rhovec[low_surv_pos] : maximum(rhovec) #missing
    # mean_rhomin = mean(rhomin_values)
    # return mean_rhomin

end