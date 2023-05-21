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

function lng(z)
    m = 20
    b = Array{Float64}(undef, m+1)
    b[1] = 1.0
    b[2] = -0.5
    for n in 2:m
        if n%2 == 0
            local sum = 0
            for k in 0:n-1
                sum += (factorial(n)/(factorial(k)*factorial(n-k)))*(b[k+1]/(n-k+1))
            end
            b[n+1] = -sum
        else
            b[n+1] = 0
        end
    end
    b = b[3:2:m+1]

    sum = 0
    for k in 1:10
        sum += b[k]*z^(1-2*k)/((2*k)*(2*k-1))
    end
    return (z-1/2)*(log(z))-z+(1/2)*log(2*π)+sum
end

function Bessel(x,z)
    sum = 0
    n=20
    for i in 0:n
        sum+= (((-1)^i)/(gamma(i+z+1)*factorial(i)))*(x/2)^(2*i)
    end
    return ((x/2)^z)*sum
end

function F11(a,b,x)
    sum = 0
    m = 150
    for n in 0:m
        sum += gamma(a+n)*exp(-big(lng(b+n)))*(big(x)^n/factorial(big(n)))
    end
    return (exp(big(lng(b)))/gamma(a))*sum
end

x1(t) = (chi_p / omega_p)*(1-exp(-omega_p*t))
y1(t) = -(c1/(2*c2))*t+(m/c2)*log(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))

u(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
v(t) = (v2_p/3)*t+(4*m/(3*c2))*log(imag(conj(r)*exp((im*eta_p/m)*(1-exp(-phi_p*t)))*exp(-im*mu_p*phi_p*t)*eta_p*exp(-phi_p*t/2)*F11(im*mu_p-im*kappa_p+1/2,1+2*im*mu_p,im*eta_p*exp(-phi_p*t))))
x2(t) = (v(t-t1)-u(t-t1))/2 + d2_x
y2(t) = (v(t-t1)+u(t-t1))/2 + d2_y

function quadDragAprox(T ; track = true)
    x = []
    y = []
    tr = 0
    cnt = 0
    n = length(T)
    for t in T
        push!(x,x1(t))
        push!(y,y1(t))
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

function instInputs(;theta = .95 , velocity = 5.0 , mass = .1 , diameter = .1)
    #input values
    global θ = theta
    global v0 = velocity
    global m = mass
    global D = diameter

    global g = 9.8
    global linC = .0016
    global quadC = .25
    global c1 = linC*D
    global c2 = quadC*D^2

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
    global k = -(π/(omega_p*sinh(π*zeta)))*(Bessel(xi_p , im*zeta)*((c2*v0*sinpi(θ/2)/m)+(c1/(2*m))-im*zeta*omega_p)+xi_p*omega_p*Bessel(xi_p,im*zeta-1))

    global d2_x = x1(t1)
    global d2_y = y1(t1)
    global v2_p = (y1(t1+h)-y1(t1-h))/(2*h)+(x1(t1+h)-x1(t1-h))/(2*h)
    global v2_m = (y1(t1+h)-y1(t1-h))/(2*h)-(x1(t1+h)-x1(t1-h))/(2*h)
    global phi_p = (2*c1+c2*v2_p)/(2*m)
    global lambda_p = v2_m + g/phi_p
    global eta_p = (3*c2*lambda_p)/(2*m*phi_p)
    global kappa_p = (3*c2*g)/(4*m*phi_p^2)
    global mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
    global r = -(1/(2*eta_p*mu_p))*(F11(im*mu_p-im*kappa_p+1/2,1+2*im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*F11(im*mu_p-im*kappa_p-1/2,1+2*im*mu_p,im*eta_p))

    
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
    println("phi_p = ",phi_p)
    println("lambda_p = ",lambda_p)
    println("eta_p = ",eta_p)
    println("kappa_p = ",kappa_p)
    println("mu_p = ",mu_p)
    println("r = ",r)
    println("\n")
end
