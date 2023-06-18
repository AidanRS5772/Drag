using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

tick()
errorData = errorTpV1(25,200,100)
tock()
println("\n\n Error Data: \n\n")
println(errorData)



