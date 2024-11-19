function breakpoint_find(massvec, rhominvec)
    # Transform to log-log space
    log_massvec = log.(massvec)
    log_rhominvec = log.(rhominvec)

    # Initialize variables to track the best breakpoint
    best_error = Inf
    best_index = nothing

    # Loop through possible breakpoints
    for i in 2:length(massvec)-2  # Avoid endpoints as breakpoints
        # Split the data
        log_massvec1, log_rhominvec1 = log_massvec[1:i], log_rhominvec[1:i]
        log_massvec2, log_rhominvec2 = log_massvec[i+1:end], log_rhominvec[i+1:end]

        # Fit linear models
        X1 = hcat(ones(length(log_massvec1)), log_massvec1)
        coeffs1 = X1 \ log_rhominvec1
        error1 = sum((X1 * coeffs1 .- log_rhominvec1).^2)

        X2 = hcat(ones(length(log_massvec2)), log_massvec2)
        coeffs2 = X2 \ log_rhominvec2
        error2 = sum((X2 * coeffs2 .- log_rhominvec2).^2)

        # Total error
        total_error = error1 + error2

        # Update the best breakpoint
        if total_error < best_error
            best_error = total_error
            best_index = i
        end
    end

    # Return the best breakpoint and its index
    return massvec[best_index], best_index
end
