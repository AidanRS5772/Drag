using PlotlyJS
using Interpolations
using Printf

include("projwDrag.jl")

function pProjP()
    sim = qDSim(track = false)
    simX = sim[1]
    simY = sim[2]
    simT = sim[3]

    ap = qDAproxP()
    apX = ap[1]
    apY = ap[2]
    apT = ap[3]

    ptb = scatter(x = [x1(t1) , x2(t2) , x3(t3) , x4(t4)] , y = [y1(t1) , y2(t2) , y3(t3) , y4(t4)] , mode = "markers" , name = "Transition Points Both")
    pty = scatter(x = [t1 , t2 , t3 , t4] , y = [y1(t1) , y2(t2) , y3(t3) , y4(t4)] , mode = "markers" , name = "Transition Points Y")
    ptx = scatter(x = [t1 , t2 , t3 , t4] , y = [x1(t1) , x2(t2) , x3(t3) , x4(t4)] , mode = "markers" , name = "Transition Points X")
    ps = scatter(x=simX , y=simY , mode = "line" , name = "Sim")
    pa = scatter(x = apX , y=apY , mode = "line" , name = "Aprox")
    psx = scatter(x = simT , y = simX , mode = "line" , name = "Sim X")
    psy = scatter(x = simT , y = simY , mode = "line" , name = "Sim Y")
    pax = scatter(x = apT , y = apX , mode = "line" , name = "Aprox X")
    pay = scatter(x = apT , y = apY, mode = "line" , name = "Aprox Y")

    return [ps , pa , ptb] , [psx , psy , pax , pay , pty , ptx]
end

function pProjV()
    sim = qDSim(track = false)
    simX = sim[1]
    simY = sim[2]
    simT = sim[3]
    dts = sim[4]

    simDX = []
    simDY = []
    simDT = simT[2:end-1]
    for n in 2:length(simT)-1
        push!(simDX , (simX[n+1]-simX[n-1])/(2*dts))
        push!(simDY , (simY[n+1]-simY[n-1])/(2*dts))
    end

    ap = qDAproxV()
    apDX = ap[1]
    apDY = ap[2]
    apDT = ap[3]

    ptb = scatter(x = [vx1(t1) , vx2(t2) , vx3(t3) , vx4(t4)] , y = [vy1(t1) , vy2(t2) , vy3(t3) , vy4(t4)] , mode = "markers" , name = "Transition Points Both")
    pty = scatter(x = [t1 , t2 , t3 , t4] , y = [vy1(t1) , vy2(t2) , vy3(t3) , vy4(t4)] , mode = "markers" , name = "Transition Points Y")
    ptx = scatter(x = [t1 , t2 , t3 , t4] , y = [vx1(t1) , vx2(t2) , vx3(t3) , vx4(t4)] , mode = "markers" , name = "Transition Points X")
    ps = scatter(x=simDX , y=simDY , mode = "line" , name = "Simulation Velocity")
    pa = scatter(x = apDX , y=apDY , mode = "line" , name = "Aproximation Velocity")
    psdx = scatter(x = simDT , y = simDX , mode = "line" , name = "Sim DX")
    psdy = scatter(x = simDT , y = simDY , mode = "line" , name = "Sim DY")
    padx = scatter(x = apDT , y = apDX , mode = "line" , name = "Aprox DX")
    pady = scatter(x = apDT , y = apDY, mode = "line" , name = "Aprox DY")

    return [ps , pa , ptb] , [psdx , psdy , padx , pady , pty , ptx]
end