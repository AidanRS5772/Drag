using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(velocity = 25.0  , theta = .42)
s1 , s2 = pProjP()
p1  = plot(s1)
p2 = plot(s2)
[p1 p2]


