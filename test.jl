using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(velocity = 50.0 , theta = .63)
tick()
s1 , s2 = pProjP()
tock()
p1 = plot(s1)
p2 = plot(s2)

[p1 p2]
