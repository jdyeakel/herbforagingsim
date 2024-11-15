function gc_sample(gains,costs,probs)
    # Checked out - this works well 10/28/2024
    
    nbins = length(gains)
    flat_probs = vec(probs);
    categorical_dist = Categorical(flat_probs)
    sample_index = rand(categorical_dist)
    row, col = divrem(sample_index - 1, nbins) #subtracting one bc divrem 0-based
    selected_gain = gains[row + 1] # add one bc julia indexing is 1-based
    selected_cost = costs[col + 1]

    return selected_gain, selected_cost
end
