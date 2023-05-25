using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs(velocity = 20.0)
preCalc()
s1 , s2 = pProjV()
p1 = plot(s1)
p2 = plot(s2)
p2
#=
T = LinRange(0,1,100)
vx2_d = []
vy2_d = []
for t in T
    push!(vx2_d , vx2(t))
    push!(vy2_d , vy2(t))
end

p1 = scatter(x = T , y = vx2_d , mode = "line" , name = "Derivitive X")
p2 = scatter(x = T , y = vy2_d , mode = "line" , name = "Derivitive Y")

plot([p1 , p2])
=#