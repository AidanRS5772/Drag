using PlotlyJS
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

for i in 2:cnt-1
    push!(dx,(xs[i+1]-xs[i-1])/(2*dt))
    push!(dy,(ys[i+1]-ys[i-1])/(2*dt))
end

q = .677269

xratio = findall(abs.(dx)/abs.(dy) .< q)
yratio = findall(abs.(dy)/abs.(dx) .< q)

t1 = yratio

px = scatter(x = t,y = xs,mode = "line" , name = "x position")
py = scatter(x = t,y = ys,mode = "line" , name = "y position")
pdx = scatter(x = t[2:cnt-1],y = dx,mode = "line" , name = "x velocity")
pdy = scatter(x = t[2:cnt-1],y = dy,mode = "line" , name = "y velocity")

position = plot([px,py])
velocity = plot([pdx,pdy])

[position velocity]