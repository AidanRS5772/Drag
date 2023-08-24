using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
include("projwDrag.jl")
include("plotProjDrag.jl")
include("dataAnalysis.jl")


n = 5
angleVals = LinRange(.03 , .97 , n)
velVals = LinRange(5 , 200 , n)

avGrid = grid(angleVals , velVals , n)

ErrorData = avAnalysis(avGrid , n)

println("abs error tf : ")
display(ErrorData[1])
println("rel error tf : ")
display(ErrorData[2])

println("abs error rng : ")
display(ErrorData[3])
println("rel error rng : ")
display(ErrorData[4])

println("abs error tp : ")
display(ErrorData[5])
println("rel error tp : ")
display(ErrorData[6])

println("abs error yp : ")
display(ErrorData[7])
println("rel error yp : ")
display(ErrorData[8])

p1 = plot(surface(z = ErrorData[1] , x = angleVals , y = velVals) , Layout(title = "Absolute Error Final Time"))
p2 = plot(surface(z = ErrorData[2] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Final Time"))
p3 = plot(surface(z = ErrorData[3] , x = angleVals , y = velVals) , Layout(title = "Absolute Error Range"))
p4 = plot(surface(z = ErrorData[4] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Range"))
p5 = plot(surface(z = ErrorData[5] , x = angleVals , y = velVals) , Layout(title = "Absolute Error Peak Time"))
p6 = plot(surface(z = ErrorData[6] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Peak Time"))
p7 = plot(surface(z = ErrorData[7] , x = angleVals , y = velVals) , Layout(title = "Absolute Error Peak Height"))
p8 = plot(surface(z = ErrorData[2] , x = angleVals , y = velVals) , Layout(title = "Reletive Error Peak Height"))

[p1 p2 p3 p4 p5 p6 p7 p8]