function calculate_histogram(data, nbins)
    # Define bin edges from minimum to maximum, with nbins intervals
    min_val = minimum(data)
    max_val = maximum(data)
    bin_edges = range(min_val, max_val, length=nbins+1)
    bin_counts = zeros(Int, nbins)

    # Bin data into counts
    for val in data
        for i = 1:nbins
            if val >= bin_edges[i] && val < bin_edges[i+1]
                bin_counts[i] += 1
                break
            end
        end
    end

    # Edge case: include data that matches the maximum value
    if data[end] == max_val
        bin_counts[end] += 1
    end

    return bin_edges, bin_counts
end