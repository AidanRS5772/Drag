using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")

function rootFind(f::Function , a , b , tol)

end

q = .677269
instInputs()
preCalc(print =  false)
T = LinRange(0,10,100)
out1 = []
for t in T
    push!(out1,ryx2(t))
end

sq = scatter(x = T , y = fill(q,100) , mode = "line")
s1 = scatter(x = T , y = out1 , mode = "line")

plot([s1 , sq])