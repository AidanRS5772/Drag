using Plots
using DataFrames
include("projwDrag.jl")

q = .67726947689135146269
dataGrain = 5
v_inp = LinRange(1,45,dataGrain)
th_inp = LinRange(.05,.95,dataGrain)
D_inp = LinRange(.01,1,dataGrain)
m_inp = LinRange(.01,1,dataGrain)

Data = [0 0 0 0 0 0]
for n1 in range(1,dataGrain)
    for n2 in range(1,dataGrain)
        for n3 in range(1,dataGrain)
            for n4 in range(1,dataGrain)
                instInputs(velocity = v_inp[n1] , theta = th_inp[n2] , diameter = D_inp[n3] , mass = m_inp[n4])
                simData = quadDragSim(dt = .0001 , track =  false)
                x = simData[1]
                y = simData[2]
                cnt = simData[4]
                dt = simData[6]

                vy = []
                vx = []
                for i in 2:cnt-1
                    push!(vx,(x[i+1]-x[i-1])/(2*dt))
                    push!(vy,(y[i+1]-y[i-1])/(2*dt))
                end

                vx = abs.(vx)
                vy = abs.(vy)
                velRatio = vy./vx

                idx = findall(velRatio .<= q)
                idx = idx.*dt

                global Data = vcat(Data, [v_inp[n1] th_inp[n2] D_inp[n3] m_inp[n4] idx[1] idx[end]])
            end
        end
    end
end

DF = DataFrame(Data, :auto)
rename!(DF,["Velocity", "Theta", "Diameter", "Mass", "Tau 1", "Tau 2"])

println(DF)