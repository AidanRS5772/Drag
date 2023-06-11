using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

tick()
out = errorRV1(10.0,100.0,100)
tock()
println(out)
