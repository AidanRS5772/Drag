using SpecialFunctions
using Printf

function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting, this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*a-(c2/m)*a*sqrt(a^2+b^2)
end

function qDSim(;dt = 10.0^-6,track = true)
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

    return x , y , time , cnt
end

function newtonsWp(f::Function , g::Function , init , tol , preP::Bool , q)
    h = 10^(-12)
    px = init

    while tol < abs(g(px))
        px -= 2*g(px)*h/(g(px+h)-g(px-h))
    end

    if preP
        x = px - tol
    else
        x = px + tol
    end
    
    while tol < abs(f(x)-q)
        x -= 2*(f(x)-q)*h/(f(x+h)-f(x-h))
        if (x > px && preP) || (x < px && !preP)
            x = 2*px - x
        end
    end
    return x
end

function newtonNp(f::Function , init , tol , q)
    h = 10^(-12)
    x = init
    while tol < abs(f(x)-q)
        x -= 2*(f(x)-q)*h/(f(x+h)-f(x-h))
    end
    return x
end

function finiteD(x , y , dt)
    dx = (x[3:end] .- x[1:end-2])/(2*dt)
    dy = (y[3:end] .- y[1:end-2])/(2*dt)
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

function GammaAll(x)
    if isinteger(x)
        x = convert(Int128,x)
        if x <= 21
            return factorial(x-1)
        elseif x <= 0
            return Inf64
        else
            return factorial(big(x-1))
        end
    elseif typeof(x) == Float64
        if x > 0
            out  = gamma(x)
            if isinf(out^2)
                return exp(big(lng(x)))
            else
                return out
            end
        else
            out = -π/(sin(π*x)*gamma(1-x))
            if out^2 != 0.0
                return out
            else
                return -π*exp(-big(lng(1-x)))/sin(π*x)
            end
        end
    elseif typeof(x) == ComplexF64
        out = gamma(x)
        if isinf(out^2) || out^2 == 0.0
            return exp(big(lng(x)))
        else
            return out
        end
    end
end

function Bessel(x,z)
    sum = 0
    n=20
    for i in 0:n
        sum+= (((-1)^i)/(GammaAll(i+z+1)*GammaAll(i+1)))*(x/2)^(2*i)
    end
    return ((x/2)^z)*sum
end

function rWhittaker(a,b,x)
    sum = 0
    m = 1200
    for n in 0:m
        sum += (GammaAll(b-a+1/2+n)/GammaAll(2*b+1+n))*(big(x)^n/GammaAll(n+1))
    end
    return convert(ComplexF64 , sqrt(x)*exp(-x/2)*(GammaAll(2*b+1)/GammaAll(b-a+1/2))*sum)
end

x1(t) = (chi_p / omega_p)*(1-exp(-omega_p*t))
y1(t) = -(c1/(2*c2))*t+(m/c2)*log(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))
vx1(t) = chi_p*exp(-omega_p*t)
vy1(t) = -(c1/(2*c2))+(m*xi_p*omega_p/(2*c2))*exp(-omega_p*t)*((imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta+1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))-(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta-1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))
rxy(t) = vx1(t)/vy1(t)

u2(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
v2(t) = (v2_p/3)*t+(4*m/(3*c2))*log(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t))))

x2(t) = (v2(t-t1)-u2(t-t1))/2 + d2_x
y2(t) = (v2(t-t1)+u2(t-t1))/2 + d2_y

x3(t) = (v3_x/2)*(t-t2)+(m/c2)*log(imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))+d3_x
y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t-t2)))-(g/omega_0)*(t-t2)+d3_y

v4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t))))
x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y

x5(t) = (chi_m/omega_m)*(1-exp(-omega_m*(t-t4))) + d5_x
y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))+d5_y

function qDAproxP(T ; track = true)
    preCalc()
    x = []
    y = []
    tr = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < t1
            push!(x,x1(t))
            push!(y,y1(t))
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

function qDAproxV(T ; track = true)
    preCalc()
    vx = []
    vy = []
    tr = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < t1
            push!(vx,vx1(t))
            push!(vy,vy1(t))
        end

        if track
            cnt += 1
            if floor(Int,100*cnt/n) > tr
                tr += 1
                @printf("%.f%% \n",100*cnt/n)
            end
        end
    end

    return vx , vy
end

function instInputs(;theta = .95 , velocity = 5.0 , mass = .1 , diameter = .1)
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
end

function preCalc(;print = true)
    q = .677269
    h = 10^(-12)

    global chi_p = v0*cospi(θ/2)
    global omega_p = (c1 + c2 * v0 * sinpi(θ/2)) / m
    global xi_p = (c2 * chi_p) / (sqrt(2) * omega_p * m)
    global zeta = sqrt(4*g*c2*m - c1^2) / (2*m*omega_p)
    global k = -(π/(omega_p*sinh(π*zeta)))*(Bessel(xi_p , im*zeta)*((c2*v0*sinpi(θ/2)/m)+(c1/(2*m))-im*zeta*omega_p)+xi_p*omega_p*Bessel(xi_p,im*zeta-1))
    
    global t1 = newtonsWp(rxy , vy1 , 0 , sqrt(h) , true , q)

    global d2_x = x1(t1)
    global d2_y = y1(t1)
    global v2_p = (y1(t1+h)-y1(t1-h))/(2*h)+(x1(t1+h)-x1(t1-h))/(2*h)
    global v2_m = (y1(t1+h)-y1(t1-h))/(2*h)-(x1(t1+h)-x1(t1-h))/(2*h)
    global phi_p = (2*c1+c2*v2_p)/(2*m)
    global lambda_p = v2_m + g/phi_p
    global eta_p = (3*c2*lambda_p)/(2*m*phi_p)
    global kappa_p = (3*c2*g)/(4*m*phi_p^2)
    global mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
    global r = -(1/(2*eta_p*mu_p))*(rWhittaker(im*kappa_p,im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p))
    #=
    global d3_x = x2(t2)
    global d3_y = y2(t2)
    global v3_x = (x2(t2+h)-x2(t2-h))/(2*h)
    global v3_y = (y2(t2+h)-y2(t2-h))/(2*h)
    global omega_0 = (c1+c2*v3_x)/m
    global gammay = v3_y+g/omega_0
    global xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
    global delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
    global epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
    global p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))

    global d4_x = x3(t3)
    global d4_y = y3(t3)
    global v4_p = (y3(t3+h)-y3(t3-h))/(2*h)+(x3(t3+h)-x3(t3-h))/(2*h)
    global v4_m = (y3(t3+h)-y3(t3-h))/(2*h)-(x3(t3+h)-x3(t3-h))/(2*h)
    global phi_m = (2*c1-c2*v4_m)/(2*m)
    global lambda_m = v4_p + g/phi_m
    global eta_m = (3*c2*lambda_m)/(2*m*phi_m)
    global kappa_m = (3*c2*g)/(4*m*phi_m^2)
    global mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
    global s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

    
    global d5_x = x4(t4)
    global d5_y = y4(t4)
    global v5_y = (y4(t4+h)-y4(t4-h))/(2*h)
    global chi_m = (x4(t4+h)-x4(t4-h))/(2*h)
    global omega_m = (c1-c2*v5_y)/m
    global xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
    global Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
    global l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
    global l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))
    =#
    if print
        println("")

        println("c1 = ",c1)
        println("c2 = ",c2)
        println("t1 = ",t1)

        println("")

        println("chi_p = ",chi_p)
        println("omega_p = ",omega_p)
        println("xi_p = ",xi_p)
        println("zeta = ",zeta)
        println("k = ",k)
        
        println("")

        println("phi_p = ",phi_p)
        println("lambda_p = ",lambda_p)
        println("eta_p = ",eta_p)
        println("kappa_p = ",kappa_p)
        println("mu_p = ",mu_p)
        println("r = ",r)

        println("")
        #=
        println("omega_0 = ", omega_0)
        println("gammay = ", gammay)
        println("xi_0 = ", xi_0)
        println("delta = ", delta)
        println("epsilon = ", epsilon)
        println("p = ",p)

        println("")

        println("phi_m = ", phi_m)
        println("lambda_m = ", lambda_m)
        println("eta_m = ", eta_m)
        println("kappa_m = ", kappa_m)
        println("mu_m = ", mu_m)
        println("s = ",s)

        println("")

        println("chi_m = ", chi_m)
        println("omega_m = ", omega_m)
        println("xi_m = ", xi_m)
        println("Omega = ", Omega)
        println("l_p = ",l_p)
        println("l_m = ",l_m)

        println("")
        =#
    end
end
