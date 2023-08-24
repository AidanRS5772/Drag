using PlotlyJS
using SpecialFunctions
using QuadGK
using TickTock
include("plotProjDrag.jl")
include("projwDrag.jl")

function newtonsRF(f::Function , df::Function , init)
    println("\nNewtons Root Finder:\n")

    tol = 1e-8
    x1 = init+2*tol
    x2 = init

    println("t = " , x2)
    while abs(x2-x1) > tol
        x1 = x2
        x2 = x1 - f(x1)/df(x1)
        println("t = " , x2)
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

    absErrorTFData = Matrix(undef , n , n)
    relErrorTFData = Matrix(undef , n , n)
    absErrorRangeData = Matrix(undef , n , n)
    relErrorRangeData = Matrix(undef , n , n)
    absErrorTPData = Matrix(undef , n , n)
    relErrorTPData = Matrix(undef , n , n)
    absErrorYPData = Matrix(undef , n , n)
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

            apxTf = newtonsRF(projPy , projVy , simTf)
            apxXr = projPx(apxTf)
            apxTp = secantRF(projVy , simTp; track = true)
            apxYp = projPy(apxTp)

            println("\napx tf = ",apxTf)
            println("apx rng = ",apxXr)
            println("apx tp = ",apxTp)
            println("apx yp = ",apxYp)

            absRangeError = abs(apxXr - simXr)
            relRangeError = abs(apxXr - simXr)/simXr
            absTfError = abs(apxTf - simTf)
            relTfError = abs(apxTf - simTf)/simTf
            absTpError = abs(apxTp - simTp)
            relTpError = abs(apxTp - simTp)/simTp
            absYpError = abs(apxYp - simYp)
            relYpError = abs(apxYp - simYp)/simYp

            println("\nabs error tf = ", absTfError)
            println("rel error tf = ", relTfError)
            println("abs error rng = " , absRangeError)
            println("rel error rng = " , relRangeError)
            println("abs error tp = ", absTpError)
            println("rel error tp = ", relTpError)
            println("abs error yp = ", absYpError)
            println("rel error yp = ", relYpError)

            absErrorTFData[i , j] = absTfError
            relErrorTFData[i , j] = relTfError
            absErrorRangeData[i , j] = absRangeError
            relErrorRangeData[i , j] = relRangeError
            absErrorTPData[i , j] = absTpError
            relErrorTPData[i , j] = relTpError
            absErrorYPData[i , j] = absYpError
            relErrorYPData[i , j] = relYpError
        end
    end

    return absErrorTFData , relErrorTFData , absErrorRangeData , relErrorRangeData , absErrorTPData , relErrorTPData , absErrorYPData , relErrorYPData
end
