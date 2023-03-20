using Plots
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs()
p = plotError(cutOff = .1)
display(plot(p[1],p[2]))



