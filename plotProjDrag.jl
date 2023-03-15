using Plots
using Interpolations
using Printf

include("projwDrag.jl")

function plotError(;dt = 10^-4,cutOff = .05)

    sim = quadDragSim()

    simX = sim[1]
    simY = sim[2]
    time = sim[3]
    cnt = sim[4]

    t = LinRange(0,time,floor(Int,time/dt))
    aprox = quadDragAprox(t)
    aproxX = convert.(Float64 , aprox[1])
    aproxY = convert.(Float64 , aprox[2])

    splineX = CubicSplineInterpolation(t, aproxX)
    splineY = CubicSplineInterpolation(t, aproxY)
    splineT = big.(collect(LinRange(0,time,cnt)))
    linAproxX = big.(collect(splineX(splineT)))
    linAproxY = big.(collect(splineY(splineT)))

    abserrorX = abs.(simX .- linAproxX)
    abserrorY = abs.(simY .- linAproxY)
    relerrorX = abserrorX./simX
    relerrorY = abserrorY./simY
    relT = collect(splineT)

    nrelerrorX = relerrorX[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    nrelerrorY = relerrorY[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    nrelT = relT[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]

    ep1 = plot(splineT,[abserrorX abserrorY],legend = false)
    ep2 = plot(nrelT,[nrelerrorX nrelerrorY],legend = false)
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

function plotProj(;dt = 10.0^-4)

    sim = quadDragSim()

    simX = sim[1]
    simY = sim[2]
    time = sim[3]

    t = LinRange(0,time,floor(Int,time/dt))
    aprox = quadDragAprox(t)
    aproxX = aprox[1]
    aproxY = aprox[2]

    p = plot(simX,simY,label = "Simulation")
    return plot(p,aproxX,aproxY,label = "Aproximation")
end