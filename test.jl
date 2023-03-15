using Plots
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(velocity = 10)
p1 , p2 = plotError()

display(plot(p1,p2))