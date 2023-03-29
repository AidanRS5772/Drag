using Plots
include("projwDrag.jl")
include("plotProjDrag.jl")


instInputs(velocity = 10)
p  = plotError()
p1 = plotProj()
display(plot(p1,p[1],p[2]))