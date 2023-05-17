using PlotlyJS
using Interpolations
using Printf

include("projwDrag.jl")

function plotError(;dt = 10^-3,cutOff = .05)

    sim = quadDragSim()

    simX = sim[1]
    simY = sim[2]
    time = sim[3]
    cnt = sim[4]

    t = LinRange(0,time,floor(Int,time/dt))
    aprox = quadDragAprox(t)
    aproxX = convert.(Float64,aprox[1])
    aproxY = convert.(Float64,aprox[2])

    
    splineX = CubicSplineInterpolation(t, aproxX)
    splineY = CubicSplineInterpolation(t, aproxY)
    splineT = collect(LinRange(0,time,cnt))

    aproxX = collect(splineX(splineT))
    aproxY = collect(splineY(splineT))

    abserrorX = abs.(simX .- aproxX)
    abserrorY = abs.(simY .- aproxY)

    relerrorX = abserrorX./simX
    relerrorY = abserrorY./simY
    relT = splineT

    nrelerrorX = relerrorX[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    nrelerrorY = relerrorY[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    nrelT = relT[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    
    println("Done")
    ep1 = plot(splineT,abserrorY)
    ep2 = plot(nrelT,nrelerrorY)
    return ep1 , ep2
end

function manyErrorPlots(N,b,cutOff)
    absErrors = Array{Plots.Plot, 1}(undef, N);
    relErrors = Array{Plots.Plot, 1}(undef, N);

    for i in 1:N
        dt = 10.0^(-i-b+1)
        @printf("\n\n\n!!New Plot!!   dt:%.6f \n\n\n",dt)
        error = plotError(dt,cutOff)
        absErrors[i] = error[1]
        relErrors[i] = error[2]
    end
    error = hcat(absErrors,relErrors)
    return plot(error...,layout = (2,N))
end

function plotProj(;dt = 10.0^-3)

    sim = quadDragSim(dt = 10^-6 , track = false)

    simX = sim[1]
    simY = sim[2]
    time = sim[3]
    cnt = sim[4]
    ts = LinRange(0,time,cnt)

    t = LinRange(0,time,floor(Int,time/dt))
    aprox = quadDragAprox(t ; track = false)
    aproxX = aprox[1]
    aproxY = aprox[2]

    ps = scatter(x=simX , y=simY , mode = "line" , name = "Simulation")
    pa = scatter(x = aproxX , y=aproxY , mode = "line" , name = "Aproximation")
    return plot([ps , pa])
end