using CSV
using DataFrames
using Statistics
using GLM
using Plots

# Load the data
file_path = "/Users/justinyeakel/Dropbox/PostDoc/2024_Anna_bse/data/Fritz_Oikos_2009_data_clean.csv"
data = CSV.read(file_path, DataFrame)

# Log-transform Body Mass (BM) and Mean Particle Size (MPS)
data[!, :logBM] = log.(data.BM)
data[!, :logMPS] = log.(data.MPS)

# Fit a linear model for the whole dataset: log(MPS) = log(a) + b*log(BM)
model_all = lm(@formula(logMPS ~ logBM), data)

# Extract intercept (log(a)) and slope (b) for the whole dataset
log_a_all = coef(model_all)[1]
b_all = coef(model_all)[2]
a_all = exp(log_a_all)

println("Whole dataset:")
println("Scaling coefficient (a): $a_all")
println("Allometric exponent (b): $b_all")

# Prepare 4x4 plot layout
plot_layout = @layout [a b; c d]

# Plot for the whole dataset
p_all = scatter(data.logBM, data.logMPS, label="Data", xlabel="log(BM)", ylabel="log(MPS)", title="Whole Dataset")
plot!(p_all, data.logBM, predict(model_all), label="Fit: log(MPS) = $log_a_all + $b_all*log(BM)", color=:red)


# Fit models for each digestive type (DT)
groups = unique(data.DT)
plots = []

for group in groups
    # Subset data for the group
    group_data = filter(row -> row.DT == group, data)
    
    # Fit a linear model: log(MPS) = log(a) + b*log(BM)
    model_group = lm(@formula(logMPS ~ logBM), group_data)

    # Extract intercept (log(a)) and slope (b) for the group
    log_a_group = coef(model_group)[1]
    b_group = coef(model_group)[2]
    a_group = exp(log_a_group)

    println("\nDigestive type: $group")
    println("Scaling coefficient (a): $a_group")
    println("Allometric exponent (b): $b_group")

    # Create a plot for the group
    p_group = scatter(group_data.logBM, group_data.logMPS, label="Data", xlabel="log(BM)", ylabel="log(MPS)", title="Digestive Type: $group")
    plot!(p_group, group_data.logBM, predict(model_group), label="Fit: log(MPS) = $log_a_group + $b_group*log(BM)", color=:red)

    push!(plots, p_group)
end

# Combine all plots in a 4x4 layout
plot(p_all, plots[1], plots[2], plots[3], layout=plot_layout, size=(1200, 800))