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

    return x , y , time , cnt , dt
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

function Whittaker(a,b,x)
    sum = 0
    n = 100
    for i in 0:n
        sum += (gamma(im*b-a+1/2+i)/gamma(2*im*b+1+i))*((im^i)*(x^i)/factorial(big(i)))
    end
    sum *= gamma(1+2*im*b)/gamma(im*b-a+1/2)
    return exp(-π*b/2)*cis(-x/2)*((1+im)/sqrt(2))*cis(b*log(x))*sqrt(x)*sum
end

x1(t) = (v0*cos(θ)/ωp)*(1-exp(-ωp*t))
y1(t) = -(c1/(2*c2))*t+(m/c2)*log(real(conj(k)*BesselIm(ξp*exp(-ωp*t),im*ζ)))

x2(t) = ((m*ω0-c1)/(2*c2))*t+(m/c2)*log(imag(conj(p)*Whittaker(im*κ,μ,ξ0*exp(-ω0*t))))
y2(t) = ((q*v0*cos(θ)/ω0)*exp(-ωp*τ1)+g/(ω0^2))*(1-exp(-ω0*(t-τ1)))-(g/ω0)*(t-τ1)+y1(τ1)

x3(t) = ((q*v0*cos(θ)/ωm)*exp(-ωp*τ1-ω0*(τ2-τ1))-g/(q*ω0*ωm)*(1-exp(-ω0*(τ2-τ1))))*(1-exp(-ωm*(t-τ2)))+x2(τ2)
y3(t) = (c1/(2*c2))*t-(m/c2)*ln(lp*BesselRe(ξm*exp(-ωm*t),Ω)+lm*BesselRe(ξm*exp(-ωm*t),-Ω))


function quadDragAprox(T)
    x = []
    y = []
    track = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < τ1
            push!(x,x1(t))
            push!(y,y1(t))
        elseif τ1 <= t < τ2
            push!(x,x2(t))
            push!(y,y2(t))
        else
            push!(x,x3(t))
            push!(y,y3(t))
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

    simData = quadDragSim(dt = .0001 , track =  false)
    x = simData[1]
    y = simData[2]
    cnt = simData[4]
    dt = simData[5]

    vy = []
    vx = []
    for i in 2:cnt-1
        push!(vx,(x[i+1]-x[i-1])/(2*dt))
        push!(vy,(y[i+1]-y[i-1])/(2*dt))
    end

    q = .677269

    vx = abs.(vx)
    vy = abs.(vy)
    velRatio = vy./vx

    idx = findall(velRatio .<= q)
    idx = idx.*dt

    tau1 = idx[1]
    tau2 = idx[end]

    global omega_plus = (c1 + c2*v0*sin(theta))/m
    global chi_plus = v0*cos(theta)
    global xi_plus = (c2*chi_plus)/(sqrt(2)*omega_plus*m)
    global zeta = (m/omega_plus)*sqrt(c2*g*m - (c1^2)/4)

    global omega_0 = (c1/m) + (c2*chi_plus/m)*exp(-omega_plus*tau1)
    global gammay = q*chi_plus*exp(-omega_plus*tau1) + g/omega_0
    global xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)*exp(omega_0*tau1)
    global kappa = (c2*g)/(sqrt(2)*m*omega_0^2)
    global mu = sqrt(2*g^2*c2^2 - c1^2*omega_0^2)/(2*m*omega_0^2)

    global chi_minus = (gammay/q)*exp(-omega_0*(tau2 - tau1)) - (g/(q*omega_0))
    global omega_minus = (c1 - c2*q*chi_minus)/m
    global xi_minus = (c2*chi_minus)/(sqrt(2)*omega_minus*m)*exp(omega_minus*tau2)
    global Omega = sqrt(c2*g*m + (c1^2)/4)/(m*omega_minus)

    #=
    println("Drag Constants:")
    println("c1 = ",c1)
    println("c2 = ",c2)
    println("tau1 = ",tau1)
    println("tau2 = ",tau2)

    println("\nPositive variables:")
    println("omega_plus = ", omega_plus)
    println("chi_plus = ", chi_plus)
    println("xi_plus = ", xi_plus)
    println("zeta = ", zeta)

    println("\nZero variables:")
    println("omega_0 = ", omega_0)
    println("gammay = ", gammay)
    println("xi_0 = ", xi_0)
    println("kappa = ", kappa)
    println("mu = ", mu)

    println("\nNegative variables:")
    println("omega_minus = ", omega_minus)
    println("chi_minus = ", chi_minus)
    println("xi_minus = ", xi_minus)
    println("Omega = ", Omega)
    =#
end
