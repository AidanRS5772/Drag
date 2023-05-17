using PlotlyJS
using Printf
include("projwDrag.jl")

q = .67726947689135146269

instInputs(theta = .621)
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

xratio = abs.(dx)./abs.(dy)
yratio = abs.(dy)./abs.(dx)

map!(x -> x > 2 ? 2 : x, xratio , xratio)

xratioq = findall(xratio .> q)
yratioq = findall(yratio .< q)

xratioq = xratioq.*dt
yratioq = yratioq.*dt

t1 = xratioq[1]
t2 = yratioq[1]
t3 = yratioq[end]
t4 = xratioq[end]

println(t1)
println(t2)
println(t3)
println(t4)

px = scatter(x = t , y = xs , mode = "line" , name = "x position")
py = scatter(x = t , y = ys , mode = "line" , name = "y position")
pdx = scatter(x = t[2:cnt-1] , y = dx , mode = "line" , name = "x velocity")
pdy = scatter(x = t[2:cnt-1] , y = dy , mode = "line" , name = "y velocity")
pxratio = scatter(x = t[2:cnt-1] , y = xratio , mode = "line" , name = "dx/dy")
pyratio = scatter(x = t[2:cnt-1] , y = yratio , mode = "line" , name = "dy/dx")
pdxs = scatter(x = t[2:cnt-1] , y = dx.^2 , mode = "line" , name = "dx^2")
pdys = scatter(x = t[2:cnt-1] , y = dy.^2 , mode = "line" , name = "dy^2")
qarray = fill(q,cnt)
qplot = scatter(x = t , y = qarray , mode = "line")

position = plot([px , py])
velocity = plot([pdx , pdy])
ratio = plot([pxratio , pyratio , qplot])
square = plot([pdxs,pdys])

ratio
