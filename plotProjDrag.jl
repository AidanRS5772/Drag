using PlotlyJS
using Printf

include("projwDrag.jl")

function plotProj(; dt=2^-5, track=true)
    sim = qDSim(track=false)
    simX = sim[1]
    simY = sim[2]
    simT = sim[3]
    simVX = sim[4]
    simVY = sim[5]

    preCalc(print = track)

    apxX = [0.0]
    apxY = [0.0]
    apxT = [0.0]
    apxVX = [v0 * cospi(θ / 2)]
    apxVY = [v0 * sinpi(θ / 2)]


    t = dt
    while true
        x, y = projP(t)
        vx , vy = projV(t)

        push!(apxX, x)
        push!(apxY, y)
        push!(apxT, t)
        push!(apxVX , vx)
        push!(apxVY , vy)

        if (y < 0)
            break
        end

        if (track)
            println("t = " , t)
        end
        t += dt
    end


    if 1 > θ > 2*acot(q)/π
        tT = [t1 , t2]
        xT = [d2x , d3x]
        yT = [d2y , d3y]
        vxT = [vx1(t1) , vx2(t2)]
        vyT = [vy1(t1) , vy2(t2)]

        if d4y > 0
            push!(tT , t3)
            push!(xT , d4x)
            push!(yT , d4y)
            push!(vxT , vx3(t3))
            push!(vyT , vy3(t3))
        end

        if d5y > 0
            push!(tT , t4)
            push!(xT , d5x)
            push!(yT , d5y)
            push!(vxT , vx4(t4))
            push!(vyT , vy4(t4))
        end
    
    elseif 2*acot(q)/π >= θ > 2*atan(q)/π
        tT = [t2]
        xT = [d3x]
        yT = [d3y]
        vxT = [vx2(t2)]
        vyT = [vy2(t2)]

        if d4y > 0
            push!(tT , t3)
            push!(xT , d4x)
            push!(yT , d4y)
            push!(vxT , vx3(t3))
            push!(vyT , vy3(t3))
        end

        if d5y > 0
            push!(tT , t4)
            push!(xT , d5x)
            push!(yT , d5y)
            push!(vxT , vx4(t4))
            push!(vyT , vy4(t4))
        end
    
    elseif 2*atan(q)/π >= θ > 0
        tT = []
        xT = []
        yT = []
        vxT = []
        vyT = []

        if d4y > 0
            push!(tT , t3)
            push!(xT , d4x)
            push!(yT , d4y)
            push!(vxT , vx3(t3))
            push!(vyT , vy3(t3))
        end

        if d5y > 0
            push!(tT , t4)
            push!(xT , d5x)
            push!(yT , d5y)
            push!(vxT , vx4(t4))
            push!(vyT , vy4(t4))
        end
    end

    sSProj = scatter(x=simX, y=simY, mode="line", name="Simulation")
    sAProj = scatter(x=apxX, y=apxY, mode="line", name="Aproximation")
    sT = scatter(x=xT, y=yT, mode="markers", name="Transition Points")

    sSX = scatter(x=simT, y=simX, mode="line", name="Sim X")
    sSY = scatter(x=simT, y=simY, mode="line", name="Sim Y")
    sAX = scatter(x=apxT, y=apxX, mode="line", name="Apx X")
    sAY = scatter(x=apxT, y=apxY, mode="line", name="Apx Y")
    sTX = scatter(x=tT, y=xT, mode="markers", name="Trans X")
    sTY = scatter(x=tT, y=yT, mode="markers", name="Trans Y")

    sSVX = scatter(x=simT , y=simVX , mode = "line" , name="Sim VX")
    sSVY = scatter(x=simT , y=simVY , mode = "line" , name="Sim VY")
    sAVX = scatter(x=apxT , y=apxVX , mode = "line" , name="Apx VX")
    sAVY = scatter(x=apxT , y=apxVY , mode = "line" , name="Apx VY")
    sTVX = scatter(x=tT, y=vxT, mode="markers", name="Trans VX")
    sTVY = scatter(x=tT, y=vyT, mode="markers", name="Trans VY")

    return [sSProj, sAProj, sT], [sSX, sAX, sSY, sAY, sTX, sTY] , [sSVX , sSVY , sAVX , sAVY , sTVX , sTVY]
end
