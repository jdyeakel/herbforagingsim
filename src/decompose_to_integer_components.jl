function decompose_to_integer_components(value::Float64)
    components = Float64[]
    remaining = value
    
    # Subtract 1 until remaining is less than 1
    while remaining >= 1
        push!(components, 1)
        remaining -= 1
    end
    
    # Add the final fractional remainder
    push!(components, remaining)
    
    return components
end