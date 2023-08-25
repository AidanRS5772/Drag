using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
include("projwDrag.jl")
include("plotProjDrag.jl")
include("dataAnalysis.jl")

instInputs(angle = .99 , velocity = 20.0)
s1, s2, s3 = plotProj()
p1 = plot(s1)
p2 = plot(s2)
p3 = plot(s3)
[p1 p2 p3]
