using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

tick()
out = errorTiTh1(.625,.95,36)
tock()
println(out)
