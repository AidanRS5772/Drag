using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(diameter = .3)

tick()
s1 , s2 = pProjP(dt = 10^(-2))
tock()
p1 = plot(s1)
p2 = plot(s2)
[p1 p2]