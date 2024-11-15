function findrhomin(rhoexpvec,probsurv)
    
    rhovec = 10 .^ rhoexpvec;

    first_nonzero_pos = findfirst(x -> x > 0, probsurv)

    # Retrieve the corresponding value from `rhovec`
    rhomin = first_nonzero_pos !== nothing ? rhovec[first_nonzero_pos] : maximum(rhovec) #missing


    return rhomin

end