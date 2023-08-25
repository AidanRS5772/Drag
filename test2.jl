using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
using DataFrames
using CSV
include("projwDrag.jl")
include("plotProjDrag.jl")
include("dataAnalysis.jl")


n = 100
angleVals = LinRange(.01 , .99 , n)
velVals = LinRange(10 , 200 , n)

a_label = ["θ = $(x)" for x in angleVals]
v0_label = ["v0 = $(x)" for x in velVals]

avGrid = grid(angleVals , velVals , n)

println(avGrid)

ErrorData = avAnalysis(avGrid , n)

println("\nReletive Error Time of Flight : \n")
tfData = DataFrame(ErrorData[1] , :auto)
rename!(tfData , a_label)
tfData = hcat(DataFrame(Velocity_Vals = v0_label), tfData)

CSV.write("Data/RE_timeFlight.csv" , tfData)
println(tfData)

println("\nReletive Error Range : \n")
rngData = DataFrame(ErrorData[2] , :auto)
rename!(rngData , a_label)
rngData = hcat(DataFrame(Velocity_Vals = v0_label), rngData)

CSV.write("Data/RE_range.csv" , rngData)
println(rngData)

println("\nReletive Error Time to Peak : \n")
tpData = DataFrame(ErrorData[3] , :auto)
rename!(tpData , a_label)
tpData = hcat(DataFrame(Velocity_Vals = v0_label), tpData)

CSV.write("Data/RE_timePeak.csv" , tpData)
println(tpData)

println("\nReletive Error Peak Height : \n")
ypData = DataFrame(ErrorData[4] , :auto)
rename!(ypData , a_label)
ypData = hcat(DataFrame(Velocity_Vals = v0_label), ypData)

CSV.write("Data/RE_heightPeak.csv" , ypData)
println(ypData)

p1 = plot(surface(z = ErrorData[1] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Final Time"))

p2 = plot(surface(z = ErrorData[2] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Range"))

p3 = plot(surface(z = ErrorData[3] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Peak Time"))

p4 = plot(surface(z = ErrorData[4] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Peak Height"))


[p1 p2 ; p3 p4]

