using PlotlyJS
using SpecialFunctions
using QuadGK
using TickTock
include("plotProjDrag.jl")
include("projwDrag.jl")

function newtonsRF(f::Function , df::Function , init ; track = true)
    if (track) println("\nNewtons Root Finder:\n") end

    tol = 1e-8
    x1 = init+2*tol
    x2 = init

    if (track) println("t = " , x2) end
    while abs(x2-x1) > tol
        x1 = x2
        x2 = x1 - f(x1)/df(x1)
        if (track) println("t = " , x2) end
    end

    return x2
end

function grid(A , B , n)
    out = Matrix(undef, n, n)

    for i in 1:n
        for j in 1:n
            out[i , j] = (A[i] , B[j])
        end
    end

    return out
end

function simErrorAnalysis(a , b)

    diffData = []

    sim1 = qDSim(;dt = 2.0^-a)
    simX1 = sim1[1]
    x1 = simX1[end]

    for n in a+1:b
        println("2^$(-n) : ")
        tick()
        sim = qDSim(;dt = 2.0^-n)
        tock()
        simX = sim[1]
        x2 = simX[end]
        
        push!(diffData , abs(x2-x1))
        x1 = x2
    end

    orderVals = a:(b-1)
    s = scatter(x = orderVals , y = diffData , mode = "markers")
    return s
end

function avAnalysis(G , n)
    time_step = 2.0^-18

    relErrorTFData = Matrix(undef , n , n)
    relErrorRangeData = Matrix(undef , n , n)
    relErrorTPData = Matrix(undef , n , n)
    relErrorYPData = Matrix(undef , n , n) 

    for i in 1:n
        for j in 1:n
            theta , v = G[i , j]
            println("\n\n\nNew Sim:\nθ = " , theta , " , v0 = " , v)

            instInputs(angle = theta , velocity = v)

            simData = qDSim(dt = time_step)

            simXr = simData[1][end]
            simTf = simData[3][end]
            simYp , idxYp = findmax(simData[2])
            simTp = simData[3][idxYp]

            println("\nsim rng = " , simXr)
            println("sim tf = " , simTf)
            println("sim yp = " , simYp)
            println("sim tp = " , simTp)


            preCalc(print = false)

            apxTf = newtonsRF(projPy , projVy , simTf ; track = false)
            apxXr = projPx(apxTf)
            apxTp = secantRF(projVy , simTp; track = false)
            apxYp = projPy(apxTp)

            println("\napx tf = ",apxTf)
            println("apx rng = ",apxXr)
            println("apx tp = ",apxTp)
            println("apx yp = ",apxYp)

            relRangeError = abs(apxXr - simXr)/simXr
            relTfError = abs(apxTf - simTf)/simTf
            relTpError = abs(apxTp - simTp)/simTp
            relYpError = abs(apxYp - simYp)/simYp

            println("\nrel error tf = ", relTfError)
            println("rel error rng = " , relRangeError)
            println("rel error tp = ", relTpError)
            println("rel error yp = ", relYpError)

            
            relErrorTFData[i , j] = relTfError
            relErrorRangeData[i , j] = relRangeError
            relErrorTPData[i , j] = relTpError
            relErrorYPData[i , j] = relYpError
        end
    end

    return relErrorTFData , relErrorRangeData , relErrorTPData , relErrorYPData
end
