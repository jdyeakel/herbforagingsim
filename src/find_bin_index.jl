function find_bin_index(value, bins)
    nbins = length(bins) - 1
    # Loop through the bins to find where the value fits.
    for i = 1:nbins
        if value >= bins[i] && value < bins[i+1]
            return i
        end
    end
    # If the value is exactly equal to the maximum, assign it to the last bin.
    return nbins
end