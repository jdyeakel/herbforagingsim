using CSV
using DataFrames
using GLM
using Plots
using Statistics

file_path = joinpath(homedir(), "Dropbox", "PostDoc", "2024_Anna_bse", "data", "OwenSmith_foragingtime.csv")

df = CSV.read(file_path, DataFrame)

df.log_foraging = log.(df.foraging_percent_day)
df.log_mass = log.(df.mass_kg)
model = lm(@formula(log_foraging ~ log_mass), df)
# Extract coefficients
intercept_log_a = coef(model)[1]  # This is log(a)
b = coef(model)[2]                # This is b

# Convert intercept to a by exponentiating
a = exp(intercept_log_a)

println("a = $a")
println("b = $b")

mass_range = range(minimum(df.mass_kg), maximum(df.mass_kg), length=100)
fit_line = a .* mass_range .^ b
Plots.scatter(df.mass_kg, df.foraging_percent_day, label="Data", xlabel="Mass (kg)", ylabel="Foraging Percent (%)", title="Foraging Percent vs. Mass", xscale=:log10, yscale=:log10)
Plots.plot!(mass_range, fit_line, label="Fit: foraging_percent = a * M^b", lw=2)


# Calculate standard errors and 95% confidence intervals
coef_table = DataFrame(coeftable(model))
intercept_se = coef_table[!,:"Std. Error"][1]
slope_se = coef_table[!,:"Std. Error"][2]

# Compute 95% confidence intervals
intercept_ci_log = (intercept_log_a - 1.96 * intercept_se, intercept_log_a + 1.96 * intercept_se)
a_ci = (exp(intercept_ci_log[1]), exp(intercept_ci_log[2]))  # Convert log(a) CI to a CI
b_ci = (b - 1.96 * slope_se, b + 1.96 * slope_se)

# Print confidence intervals
println("95% CI for a: $a_ci")
println("95% CI for b: $b_ci")