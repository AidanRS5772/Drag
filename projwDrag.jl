using SpecialFunctions
using Printf

function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting, this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*a-(c2/m)*a*sqrt(a^2+b^2)
end

function quadDragSim(;dt = 10.0^-6,track = true)
    if track
        println("\nSimulation:\n")
    end

    time = 0
    cnt = 0
    tp = 0
    y = []
    x = []

    vy=v0*sinpi(θ/2)
    vx=v0*cospi(θ/2)

    ny = 0
    nx = 0

    while ny > -eps(Float64) 

        time += dt
        cnt += 1

        if floor(Int,time*1000)>tp
            tp += 1
            if track
                println("$(tp) ms")
            end
        end

        push!(x,nx)
        push!(y,ny)

        xk1 = dragEq(vx,vy)
        yk1 = dragEq(vy,vx)-g

        xk2 = dragEq(vx+dt*xk1/2,vy+dt*yk1/2)
        yk2 = dragEq(vy+dt*yk1/2 , vx+dt*xk1/2)-g

        xk3 = dragEq(vx+dt*xk2/2 , vy+dt*yk2/2)
        yk3 = dragEq(vy+dt*yk2/2 , vx+dt*xk2/2)-g

        xk4 = dragEq(vx+dt*xk3 , vy+dt*yk3)
        yk4 = dragEq(vy+dt*yk3 , vx+dt*xk3)-g

        vx += (xk1+2*xk2+2*xk3+xk4)*dt/6
        vy += (yk1+2*yk2+2*yk3+yk4)*dt/6
        
        nx = x[end]+dt*vx
        ny = y[end]+dt*vy
    end

    ymaxi = findmax(y)[2]

    return x , y , time , cnt , ymaxi*time/cnt , dt
end

function BesselIm(x,z)
    sum = 0
    n=20
    for i in 0:n
        sum+= (((-1)^i)/(factorial(i)*gamma(i+z+1)))*(x/2)^(2*i)
    end
    return ((x/2)^z)*sum
end

function BesselRe(x,z)
    sum = 0
    n=20
    for i in 0:n
        if i+z+1 >= 0
            g = gamma(big(i+z+1))
        else
            g = (-1)^(i)*π/(sinpi(z)*gamma(-big(z+i)))
        end
        sum+= (((-1)^i)/(factorial(i)*g))*(x/2)^(2*i)
    end
    return ((big(x)/2)^z)*sum
end

function nx1(t)
    sum = 0
    for n in 0:20
        term1 = (ξm/2)^(im*ζ)*conj(k)/((2*n+im*ζ+σ/ωm)*gamma(im*ζ+n+1))
        term2 = imag(term1)-exp(-(σ+2*n*ωm)*t)*imag(cis(-ζ*ωm*t)*term1)
        sum += term2*((-1)^n)*(ξm^(2*n))/(factorial(n)*4^n)
    end
     
    mult = v0*sinpi(θ/2)/(imag(conj(k)*BesselIm(ξm,im*ζ))*ωm)
    return mult*sum
end

x1(t) = (v0*cospi(θ/2)/ωm)*(1-exp(-ωm*t))
y1(t) = -(c1/(2*c2))*t+(m/c2)*log(-imag(conj(k)*BesselIm(ξm*exp(-ωm*t),im*ζ)))+(m/c2)*log(Ξ1)

x2(t) = -(v0*cospi(θ/2)/ωp)*exp((ωp-ωm)*τ)*exp(-ωp*t)+v0*cospi(θ/2)*((1/ωp-1/ωm)*exp(-ωm*τ)+1/ωm)
y2(t) = (c1/(2*c2))*t-(m/c2)*log(-l1*BesselRe(ξp*exp(-ωp*t),Ω)+l2*BesselRe(ξp*exp(-ωp*t),-Ω))

function quadDragAprox(T)
    x = []
    y = []
    track = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < τ
            push!(x,x1(t))
            push!(y,y1(t))
        else
            push!(x,x2(t))
            push!(y,y2(t))
        end

        cnt += 1
        if floor(Int,100*cnt/n) > track
            track += 1
            @printf("%.f%% \n",100*cnt/n)
        end
    
    end

    return x , y
end

function instInputs(;theta = .95 , velocity = 1.0 , mass = .1 , diameter = .05)
    #input values
    global θ = theta
    global v0 = velocity
    global m = mass
    global D = diameter

    global g = 9.8
    global linC = .00016
    global quadC = .25
    global c1 = linC*D
    global c2 = quadC*D^2

    simData = quadDragSim(dt = .0001,track = false)
    tp = simData[5]
    ep = simData[3]*10^-3
    cnt = simData[4]
    dt = simData[6]
    global τ = tp + ep

    global am = v0*sinpi(θ/2)
    global ap = (simData[2][floor(Int,τ/dt)+1]-simData[2][floor(Int,τ/dt)-1])/(2*dt)

    global ζ = sqrt(c2*g*m-(c1/2)^2)/(c1+c2*am)
    global ξm = c2*v0*cospi(θ/2)/(sqrt(2)*(c1+c2*am))
    global ωm = (c1+c2*am)/m
    global σ = 3*c1/(2*m)
    global k =  BesselIm(ξm,im*ζ)*((2*c2*v0*sinpi(θ/2)+c1)/(m*ξm*ωm))+BesselIm(ξm,im*ζ-1)-BesselIm(ξm,im*ζ+1)
    global Ξ1 = π*ξm/(2*sinh(π*ζ))

    global Ω = sqrt(c2*g*m+(c1/2)^2)/(c1-c2*ap)
    global ωp = (c1-c2*ap)/m
    global ξp = c2*v0*cospi(θ/2)*exp((ωp-ωm)*τ)/(sqrt(2)*(c1-c2*ap))
    global Λ = (2*c1*exp(ωm*τ)/(m*ωp*ξp))+exp((ωp-ωm)*τ)*(ωm*ξm/(ωp*ξp))*(imag(conj(k)*(BesselIm(ξm*exp(-ωm*τ),im*ζ-1)-BesselIm(ξm*exp(-ωm*τ),im*ζ+1)))/imag(conj(k)*(BesselIm(ξm*exp(-ωm*τ),im*ζ))))
    global Ξ2 = π*ξp*exp((c1/m-ωp)*τ)/(4*sinpi(Ω)*Ξ1*imag(conj(k)*BesselIm(ξm*exp(-ωm*τ),im*ζ)))
    global l1 = (BesselRe(ξp*exp(-ωp*τ),-Ω)*Λ+BesselRe(ξp*exp(-ωp*τ),-Ω-1)-BesselRe(ξp*exp(-ωp*τ),-Ω+1))*Ξ2
    global l2 = (BesselRe(ξp*exp(-ωp*τ),Ω)*Λ+BesselRe(ξp*exp(-ωp*τ),Ω-1)-BesselRe(ξp*exp(-ωp*τ),Ω+1))*Ξ2
end
