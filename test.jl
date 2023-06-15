using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(velocity = 60.0  , theta = .9)
s1 , s2 = pProjV()
p1  = plot(s1)
p2 = plot(s2)
[p1 p2]


