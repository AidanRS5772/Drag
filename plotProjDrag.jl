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

function pProjP(;dt = 10.0^-3)
    sim = qDSim(track = false)

    simX = sim[1]
    simY = sim[2]
    time = sim[3]
    cnt = sim[4]

    ta = LinRange(0,time,floor(Int,time/dt))
    ts = LinRange(0,time,cnt)
    aprox = qDAproxP(ta ; track = true)
    aproxX = aprox[1]
    aproxY = aprox[2]

    ptb = scatter(x = [x1(t1) , x2(t2) , x3(t3) , x4(t4)] , y = [y1(t1) , y2(t2) , y3(t3) , y4(t4)] , mode = "markers" , name = "Transition Points Both")
    pty = scatter(x = [t1 , t2 , t3 , t4] , y = [y1(t1) , y2(t2) , y3(t3) , y4(t4)] , mode = "markers" , name = "Transition Points Y")
    ptx = scatter(x = [t1 , t2 , t3 , t4] , y = [x1(t1) , x2(t2) , x3(t3) , x4(t4)] , mode = "markers" , name = "Transition Points X")
    ps = scatter(x=simX , y=simY , mode = "line" , name = "Sim")
    pa = scatter(x = aproxX , y=aproxY , mode = "line" , name = "Aprox")
    psx = scatter(x = ts , y = simX , mode = "line" , name = "Sim X")
    psy = scatter(x = ts , y = simY , mode = "line" , name = "Sim Y")
    pax = scatter(x = ta , y = aproxX , mode = "line" , name = "Aprox X")
    pay = scatter(x = ta , y = aproxY, mode = "line" , name = "Aprox Y")

    return [ps , pa , ptb] , [psx , psy , pax , pay , pty , ptx]
end

function pProjV(;dt = 10.0^-3)
    dts = 10^(-5)
    sim = qDSim(dt = dts , track = false)
    simX = sim[1]
    simY = sim[2]
    time = sim[3]
    cnt = sim[4]

    simDX = []
    simDY = []
    for n in 2:cnt-1
        push!(simDX , (simX[n+1]-simX[n-1])/(2*dts))
        push!(simDY , (simY[n+1]-simY[n-1])/(2*dts))
    end

    ta = LinRange(0,time,floor(Int,time/dt))
    ts = LinRange(0,time,cnt)
    aproxV = qDAproxV(ta ; track = true)
    aproxDX = aproxV[1]
    aproxDY = aproxV[2]

    ptb = scatter(x = [vx1(t1) , vx2(t2) , vx3(t3) , vx4(t4)] , y = [vy1(t1) , vy2(t2) , vy3(t3) , vy4(t4)] , mode = "markers" , name = "Transition Points Both")
    pty = scatter(x = [t1 , t2 , t3 , t4] , y = [vy1(t1) , vy2(t2) , vy3(t3) , vy4(t4)] , mode = "markers" , name = "Transition Points Y")
    ptx = scatter(x = [t1 , t2 , t3 , t4] , y = [vx1(t1) , vx2(t2) , vx3(t3) , vx4(t4)] , mode = "markers" , name = "Transition Points X")
    ps = scatter(x=simDX , y=simDY , mode = "line" , name = "Simulation Velocity")
    pa = scatter(x = aproxDX , y=aproxDY , mode = "line" , name = "Aproximation Velocity")
    psdx = scatter(x = ts , y = simDX , mode = "line" , name = "Sim DX")
    psdy = scatter(x = ts , y = simDY , mode = "line" , name = "Sim DY")
    padx = scatter(x = ta , y = aproxDX , mode = "line" , name = "Aprox DX")
    pady = scatter(x = ta , y = aproxDY, mode = "line" , name = "Aprox DY")

    return [ps , pa , ptb] , [psdx , psdy , padx , pady , pty , ptx]
end