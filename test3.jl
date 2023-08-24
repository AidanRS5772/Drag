using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
using Random
include("projwDrag.jl")
include("plotProjDrag.jl")
include("dataAnalysis.jl")

instInputs()
simData = qDSim()
simYp , idxYp = findmax(simData[2])
simTp = simData[3][idxYp]

println(simYp)
println(simTp)
