using Plots
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(velocity = 10 , diameter = .5)
display(plotProj())