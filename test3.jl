using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
using Random
using DataFrames
using CSV
include("projwDrag.jl")
include("plotProjDrag.jl")
include("dataAnalysis.jl")

n = 5
mat = rand(5, 5)
df = DataFrame(mat , :auto)

col_label = ["col = $x" for x in 1:n]
row_label = ["row = $x" for x in 1:n]
rename!(df , col_label)
df = hcat(DataFrame(RowNames = row_label), df)
CSV.write("Data/test.csv", df)
println(df)
