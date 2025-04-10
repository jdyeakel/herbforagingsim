using CSV
using DataFrames
using GLM
using Plots
using Statistics

file_path = joinpath(homedir(), "Dropbox", "PostDoc", "2024_herbforaging","herbforagingsim", "data", "pawar_reactiondistance.csv")

#Log10 values
df = CSV.read(file_path, DataFrame)
data = Array(df)

massvec = collect(0.01:1:1000)
rxnd = 5 .* massvec .^(1/3)


scatter(data[:,1],data[:,2])
plot!(log10.(massvec),log10.(rxnd))