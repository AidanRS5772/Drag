using PlotlyJS
using Interpolations
using Printf

include("projwDrag.jl")

function plotProj(; dt = 2^-6 , track = true)
    sim = qDSim(track = false)
    simX = sim[1]
    simY = sim[2]
    simT = sim[3]
    simVX = sim[4]
    simVY = sim[5]

    preCalc()
    
    apxX = [0.0]
    apxY = [0.0]
    apxT = [0.0]
    apxVX = [v0*cospi(θ/2)]
    apxVY = [v0*sinpi(θ/2)]


    t = dt
    while true
        x , y = qDAproxP(t)
        vx , vy = qDAproxV(t)

        push!(apxX , x)
        push!(apxY , y)
        push!(apxVX , vx)
        push!(apxVY , vy)
        push!(apxT , t)

        if (y < 0) break end

        if (track) println(t) end
        t += dt 
    end

    tT = []
    tX = []
    tY = []
    tVX = []
    tVY = []

    function addT1()
        push!(tT , t1)
        push!(tX , x1(t1))
        push!(tY , y1(t1))
        push!(tVX , vx1(t1))
        push!(tVY , vy1(t1))
    end

    function addT2()
        push!(tT , t2)
        push!(tX , x2(t2))
        push!(tY , y2(t2))
        push!(tVX , vx2(t2))
        push!(tVY , vy2(t2))
    end

    function addT3()
        push!(tT , t3)
        push!(tX , x3(t3))
        push!(tY , y3(t3))
        push!(tVX , vx3(t3))
        push!(tVY , vy3(t3))
    end

    function addT4()
        push!(tT , t4)
        push!(tX , x4(t4))
        push!(tY , y4(t4))
        push!(tVX , vx4(t4))
        push!(tVY , vy4(t4))
    end 

    if 1 > θ >= acot(q)*2/π
        if apxT[end] > t1
            addT1()
        end
        if apxT[end] > t2
            addT2()
        end
        if apxT[end] > t3
            addT3()
        end
        if apxT[end] > t4
            addT4()
        end
    elseif acot(q)*2/π >= θ > atan(q)*2/π
        if apxT[end] > t2
            addT2()
        end
        if apxT[end] > t3
            addT3()
        end
        if apxT[end] > t4
            addT4()
        end
    elseif atan(q)*2/π >= θ > 0
        if apxT[end] > t3
            addT3()
        end
        if apxT[end] > t4
            addT4()
        end
    end

    sSProj = scatter(x = simX , y = simY , mode = "line" , name = "Simulation")
    sAProj = scatter(x = apxX , y = apxY , mode = "line" , name = "Aproximation")
    sT = scatter(x = tX , y = tY , mode = "markers" , name = "Transition Points")
    
    sSX = scatter(x = simT , y = simX , mode = "line" , name = "Sim X")
    sSY = scatter(x = simT , y = simY , mode = "line" , name = "Sim Y")
    sAX = scatter(x = apxT , y = apxX , mode = "line" , name = "Apx X")
    sAY = scatter(x = apxT , y = apxY , mode = "line" , name = "Apx Y")
    sTX = scatter(x = tT , y = tX , mode = "markers" , name = "Trans X")
    sTY = scatter(x = tT , y = tY , mode = "markers" , name = "Trans Y")

    sSVX = scatter(x = simT , y = simVX , mode = "line" , name = "Sim VX")
    sSVY = scatter(x = simT , y = simVY , mode = "line" , name = "Sim VY")
    sAVX = scatter(x = apxT , y = apxVX , mode = "line" , name = "Apx VX")
    sAVY = scatter(x = apxT , y = apxVY , mode = "line" , name = "Apx VY")
    sTVX = scatter(x = tT , y = tVX , mode = "markers" , name = "Trans VX")
    sTVY = scatter(x = tT , y = tVY , mode = "markers" , name = "Trans VY")

    return [sSProj , sAProj , sT] , [sSX , sAX , sSY , sAY , sTX , sTY] , [sSVX , sAVX , sSVY , sAVY , sTVX , sTVY]
end