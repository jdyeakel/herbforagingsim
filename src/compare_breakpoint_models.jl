function compare_breakpoint_models(massvec, rhominvec, break_index)
    log_massvec = log.(massvec)
    log_rhominvec = log.(rhominvec)

    # Single model fit
    X_all = hcat(ones(length(log_massvec)), log_massvec)
    coeffs_all = X_all \ log_rhominvec
    error_all = sum((X_all * coeffs_all .- log_rhominvec).^2)

    # Two-piece model fit
    log_massvec1, log_rhominvec1 = log_massvec[1:break_index], log_rhominvec[1:break_index]
    log_massvec2, log_rhominvec2 = log_massvec[break_index+1:end], log_rhominvec[break_index+1:end]

    X1 = hcat(ones(length(log_massvec1)), log_massvec1)
    coeffs1 = X1 \ log_rhominvec1
    error1 = sum((X1 * coeffs1 .- log_rhominvec1).^2)

    X2 = hcat(ones(length(log_massvec2)), log_massvec2)
    coeffs2 = X2 \ log_rhominvec2
    error2 = sum((X2 * coeffs2 .- log_rhominvec2).^2)

    total_error_two_piece = error1 + error2

    p1=2;

    p2=4;   
    
    n = length(massvec);
    
    F_value = ((error_all - total_error_two_piece) / (p2 - p1)) / (total_error_two_piece / (n - p2))

    p_value = 1 - cdf(FDist(p2 - p1, n - p2), F_value)

    return error_all, total_error_two_piece, F_value, p_value
end
