using PlotlyJS
using Printf
include("projwDrag.jl")

q = .67726947689135146269

instInputs(velocity = 45 , diameter = .1)
simdata = quadDragSim(dt = .0001 , track = false)
xs = simdata[1]
ys = simdata[2]
time = simdata[3]
cnt = simdata[4]
dt = simdata[5]

t = LinRange(0,time,cnt)

dx = []
dy = []

global track = 1
for i in 2:cnt-1
    push!(dx,(xs[i+1]-xs[i-1])/(2*dt))
    push!(dy,(ys[i+1]-ys[i-1])/(2*dt))
end

q = .677269

xratio = findall(abs.(dx)./abs.(dy) .< q)
yratio = findall(abs.(dy)./abs.(dx) .< q)

xratio = xratio.*dt
yratio = yratio.*dt

t1 = yratio[1]
t2 = xratio[1]
t3 = xratio[end]
t4 = yratio[end]

println(t1)
println(t2)
println(t3)
println(t4)

px = scatter(x = t,y = xs,mode = "line" , name = "x position")
py = scatter(x = t,y = ys,mode = "line" , name = "y position")
pdx = scatter(x = t[2:cnt-1],y = dx,mode = "line" , name = "x velocity")
pdy = scatter(x = t[2:cnt-1],y = dy,mode = "line" , name = "y velocity")

position = plot([px,py])
velocity = plot([pdx,pdy])

[position velocity]