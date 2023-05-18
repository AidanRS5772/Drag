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
        sum+= (((-1)^i)/((factorial(i))*gamma(i+z+1)))*(x/2)^(2*i)
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
    n = 20
    for i in 0:n
        sum += (gamma(b-a+1/2+i)/gamma(2*b+1+i))*(x^i/factorial(i))
    end
    sum *= gamma(1+2*b)/gamma(b-a+1/2)
    return exp(-x/2)*(x^(b+1/2))*sum
end

x1(t) = (chi_p / omega_p)*(1-exp(-omega_p*t))
y1(t) = -(c1/(2*c2))*t+(m/c2)*log(imag(conj(k)*BesselIm(xi_p*exp(-omega_p*t),im*zeta)))

x2(t) = (2*m/(3*c2))*log(imag(conj(r)*Whittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*(t-t1)))))-(lambda_p/(2*phi_p))*(1-exp(-phi_p*(t-t1)))+(d_p0/6+g/(2*phi_p))*(t-t1) + d_p1
y2(t) = (2*m/(3*c2))*log(imag(conj(r)*Whittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*(t-t1)))))+(lambda_p/(2*phi_p))*(1-exp(-phi_p*(t-t1)))+(d_p0/6-g/(2*phi_p))*(t-t1) + d_p2

function quadDragAprox(T ; track = true)
    x = []
    y = []
    tr = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < t1
            push!(x,x1(t))
            push!(y,y1(t))
        else
            push!(x,x2(t))
            push!(y,y2(t))
        end

        if track
            cnt += 1
            if floor(Int,100*cnt/n) > tr
                tr += 1
                @printf("%.f%% \n",100*cnt/n)
            end
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
    global c1 = 1
    global c2 = 1

    simdata = quadDragSim(dt = .0001 , track = false)
    xs = simdata[1]
    ys = simdata[2]
    time = simdata[3]
    cnt = simdata[4]
    dt = simdata[5]

    t = LinRange(0,time,cnt)

    dx = []
    dy = []

    for i in 2:cnt-1
        push!(dx,(xs[i+1]-xs[i-1])/(2*dt))
        push!(dy,(ys[i+1]-ys[i-1])/(2*dt))
    end

    global q = .677269

    xratio = abs.(dx)./abs.(dy)
    yratio = abs.(dy)./abs.(dx)

    xratioq = findall(xratio .> q)
    yratioq = findall(yratio .< q)

    xratioq = xratioq.*dt
    yratioq = yratioq.*dt

    h = 10^(-8)

    global t1 = xratioq[1]
    global t2 = yratioq[1]
    global t3 = yratioq[end]
    global t4 = xratioq[end]


    global chi_p = v0*cospi(θ/2)
    global omega_p = (c1 + c2 * v0 * sinpi(θ/2)) / m
    global xi_p = (c2 * chi_p) / (sqrt(2) * omega_p * m)
    global zeta = sqrt(4*g*c2*m - c1^2) / (2*m*omega_p)
    global k = -(π/(omega_p*sinh(π*zeta)))*(BesselIm(xi_p , im*zeta)*((c2*v0*sinpi(θ/2)/m)+(c1/(2*m))-im*zeta*omega_p)+xi_p*omega_p*BesselIm(xi_p,im*zeta-1))

    
    println("c1 = ",c1)
    println("c2 = ",c2)
    println("t1 = ",t1)
    println("t2 = ",t2)
    println("t3 = ",t3)
    println("t4 = ",t4)
    println("chi_p = ",chi_p)
    println("omega_p = ",omega_p)
    println("xi_p = ",xi_p)
    println("zeta = ",zeta)
    println("k = ",k)
    

end
