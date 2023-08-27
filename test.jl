using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
include("projwDrag.jl")
include("plotProjDrag.jl")
include("dataAnalysis.jl")

tick()
instInputs(angle = .98 , velocity = 10)
s1 , s2 , s3 = plotProj()
tock()

[plot(s1) plot(s2) plot(s3)]