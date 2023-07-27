using PlotlyJS
using SpecialFunctions
using QuadGK
using TickTock
include("plotProjDrag.jl")
include("projwDrag.jl")

instInputs()
s1, s2, s3 = plotProj()
p1 = plot(s1)
p2 = plot(s2)
p3 = plot(s3)
[p1 p2 p3]
