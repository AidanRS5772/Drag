using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

tick()
out = errorRV3(10.0,200.0,100)
tock()
println(out)
