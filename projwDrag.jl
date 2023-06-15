using SpecialFunctions
using QuadGK
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

    cnt = 0
    tr = 0

    y = []
    x = []
    t = []

    vy=v0*sinpi(θ/2)
    vx=v0*cospi(θ/2)

    ny = 0
    nx = 0

    while ny > -eps(Float64) 

        if floor(Int,cnt*dt*1000)>tr
            tr += 1
            if track
                println("$(tr) ms")
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

        cnt += 1
    end

    time = (cnt-1)*dt
    t = 0:dt:time

    return x , y , t , dt
end

function forwardEval(f::Function , a , tol)
    fa = f(a)
    a += tol
    fb = f(a)
    while fa*fb > 0
        a += tol
        fb = f(a)
    end
    
    h = 1e-6
    p1 = a - tol
    p2 = a
    c = 0
    while abs(p1-p2) > h
        c = (p1+p2)/2
        if f(p1)*f(c) < 0 
            p2 = c
        else
            p1 = c
        end
    end

    return c
end

function numInt(f::Function , t1 , t2)
    out , ~ = quadgk( f , t1 , t2 , rtol=1e-8)
    return out
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
    n=0
    while true
        val = (((-1)^n)/(GammaAll(n+z+1)*GammaAll(n+1)))*(x/2)^(2*n )
        sum += val
        if abs(val/sum) < 10.0^(-30)
            break
        end
        n += 1
    end
    return ((x/2)^z)*sum
end

function rWhittaker(a,b,x)
    sum = 0
    n = 0
    while true
        val = (GammaAll(b-a+1/2+n)/GammaAll(2*b+1+n))*(big(x)^n/GammaAll(n+1))
        sum += val
        if abs(val/sum) < 10.0^(-30)
            break
        end
        n += 1
    end
    return convert(ComplexF64 , sqrt(x)*exp(-x/2)*(GammaAll(2*b+1)/GammaAll(b-a+1/2))*sum)
end

function qDAproxP(;dt = 10.0^(-2) , track = true)
    if θ > .621
        x , y , T , xt , yt , Tt = regime1(dt , true , track , true)
    elseif .621 >= θ >  .3789
        x , y , T , xt , yt , Tt = regime2(dt , true , track , true)
    elseif .3789 >= θ
        x , y , T , xt , yt , Tt = regime3(dt , true , track , true)
    end

    return x , y , T , xt , yt , Tt
end

function qDAproxV(;dt = 10.0^(-2) , track = true)
    if θ > .621
        x , y , T , xt , yt , Tt = regime1(dt , true , track , false)
    elseif .621 >= θ >  .3789
        x , y , T , xt , yt , Tt = regime2(dt , true , track , false)
    elseif .3789 >= θ
        x , y , T , xt , yt , Tt = regime3(dt , true , track , false)
    end

    return x , y , T , xt , yt , Tt
end

function instInputs(;theta = .95 , velocity = 10.0 , mass = .1 , diameter = .1)
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

function regime1(dt , print::Bool , track::Bool , PorV::Bool)
    #To do Position input true and velocity input false for "PorV"

    vx1(t) = chi_p*exp(-c1*t/(2*c2))/(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))
    vy1(t) = -(c1/(2*c2))+(m*xi_p*omega_p/(2*c2))*exp(-omega_p*t)*((imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta+1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))-(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta-1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))
    x1(t) = numInt(vx1 , 0, t)
    y1(t) = -(c1/(2*c2))*t+(m/c2)*log(abs(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))
    rxy1(t) = vx1(t)/vy1(t)-q

    temp_u2(t) = exp((c1/m+c2*v2_p/(6*m))*t)*((abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))^(2/3))
    du2(t) = v2_m/temp_u2(t)-(g/temp_u2(t))*numInt(temp_u2 , 0 , t)
    u2(t) = numInt(du2 , 0 , t)
    #du2(t) = lambda_p*exp(-phi_p*t)-g/phi_p
    dv2(t) = (v2_p/3)+(2*m*phi_p/(3*c2))*((2*kappa_p-eta_p*exp(-phi_p*t))*real(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))-imag(conj(r)*(1+2*im*(kappa_p+mu_p))*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p*exp(-phi_p*t))))/imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))
    #u2(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
    v2(t) = (v2_p/3)*t+(4*m/(3*c2))*log(abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))
    x2(t) = (v2(t-t1)-u2(t-t1))/2 + d2_x
    y2(t) = (v2(t-t1)+u2(t-t1))/2 + d2_y
    vx2(t) = (dv2(t-t1)-du2(t-t1))/2
    vy2(t) = (dv2(t-t1)+du2(t-t1))/2
    ryx2(t) = vy2(t)/vx2(t)-q

    temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*t)*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*t)))
    temp_vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
    vy3(t) = v3_y/temp_y3(t-t2)-(g/temp_y3(t-t2))*numInt(temp_y3 , 0 , t-t2)
    y3(t) = numInt(temp_vy3 , 0 , t-t2)+d3_y
    vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t-t2)))*real(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))
    #vy3(t) = gammay*exp(-omega_0*(t-t2))-g/omega_0
    x3(t) = (v3_x/2)*(t-t2)+(m/c2)*log(abs(imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))))+d3_x
    #y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t-t2)))-(g/omega_0)*(t-t2)+d3_y
    ryx3(t) = -vy3(t)/vx3(t)-q

    temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*((abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3))
    dv4(t) = v4_p/temp_v4(t)-(g/temp_v4(t))*numInt(temp_v4 , 0 , t)
    v4(t) = numInt(dv4 , 0 , t)
    #dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
    du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
    #v4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
    u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))
    x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
    y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
    vx4(t) = (dv4(t-t3)-du4(t-t3))/2
    vy4(t) = (dv4(t-t3)+du4(t-t3))/2
    rxy4(t) = -vx4(t)/vy4(t)-q

    vx5(t) = chi_m*exp(-c1*(t-t4)/(2*c2))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
    vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
    x5(t) = numInt(vx5 , t4 , t) + d5_x
    y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(abs(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega)))+d5_y

    q = .677269
    
    chi_p = v0*cospi(θ/2)
    omega_p = (c1 + c2 * v0 * sinpi(θ/2)) / m
    xi_p = (c2 * chi_p) / (sqrt(2) * omega_p * m)
    zeta = sqrt(4*g*c2*m - c1^2) / (2*m*omega_p)
    k = -(π/(omega_p*sinh(π*zeta)))*(Bessel(xi_p , im*zeta)*((c2*v0*sinpi(θ/2)/m)+(c1/(2*m))-im*zeta*omega_p)+xi_p*omega_p*Bessel(xi_p,im*zeta-1))
    
    t1 = forwardEval(rxy1 , 0 , .05)

    d2_x = x1(t1)
    d2_y = y1(t1)
    v2_p = vy1(t1)+vx1(t1)
    v2_m = vy1(t1)-vx1(t1)
    phi_p = (2*c1+c2*v2_p)/(2*m)
    lambda_p = v2_m + g/phi_p
    eta_p = (3*c2*lambda_p)/(2*m*phi_p)
    kappa_p = (3*c2*g)/(4*m*phi_p^2)
    mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
    r = -(1/(2*eta_p*mu_p))*(rWhittaker(im*kappa_p,im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p))
    
    t2 = forwardEval(ryx2 , t1 , .05)
    
    d3_x = x2(t2)
    d3_y = y2(t2)
    v3_x = vx2(t2)
    v3_y = vy2(t2)
    omega_0 = (c1+c2*v3_x)/m
    gammay = v3_y+g/omega_0
    xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
    delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
    epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
    p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))
    
    t3 = forwardEval(ryx3 , t2 , .05)
    
    d4_x = x3(t3)
    d4_y = y3(t3)
    v4_p = vy3(t3)+vx3(t3)
    v4_m = vy3(t3)-vx3(t3)
    phi_m = (2*c1-c2*v4_m)/(2*m)
    lambda_m = v4_p + g/phi_m
    eta_m = (3*c2*lambda_m)/(2*m*phi_m)
    kappa_m = (3*c2*g)/(4*m*phi_m^2)
    mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
    s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

    t4 = forwardEval(rxy4 , t3 , .05)
    
    d5_x = x4(t4)
    d5_y = y4(t4)
    v5_y = vy4(t4)
    chi_m = vx4(t4)
    omega_m = (c1-c2*v5_y)/m
    xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
    Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
    l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
    l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

    if print
        println("")

        println("c1 = ",c1)
        println("c2 = ",c2)
        println("t1 = ",t1)
        println("t2 = ",t2)
        println("t3 = ",t3)
        println("t4 = ",t4)

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
    end

    if PorV
        x = []
        y = []
        T = []
        tr = 0
        cnt = 0
        while true

            t = cnt*dt
            push!(T,t)
            
            if t <= t1
                push!(x,x1(t))
                push!(y,y1(t))
            elseif t1 < t <= t2
                push!(x,x2(t))
                push!(y,y2(t))
            elseif t2 < t <= t3
                push!(x,x3(t))
                push!(y,y3(t))
            elseif t3 < t <= t4
                push!(x,x4(t))
                push!(y,y4(t))
            elseif t4 < t
                push!(x,x5(t))
                push!(y,y5(t))
                if y5(t) < 0
                    break
                end
            end

            
            if floor(Int,cnt*dt*1000)>tr
                tr += 1
                if track
                    println("$(tr) ms")
                end
            end

            cnt += 1
        
        end

        return x , y , T , [x1(t1) , x2(t2) , x3(t3) , x4(t4)] , [y1(t1) , y2(t2) , y3(t3) , y4(t4)] , [t1 , t2 , t3 , t4] 
    else
        vx = []
        vy = []
        T = []
        tr = 0
        cnt = 0
        while true

            t = cnt*dt
            push!(T,t)

            if t <= t1
                push!(vx,vx1(t))
                push!(vy,vy1(t))
            elseif t1 < t <= t2
                push!(vx,vx2(t))
                push!(vy,vy2(t))
            elseif t2 < t <= t3
                push!(vx,vx3(t))
                push!(vy,vy3(t))
            elseif t3 < t <= t4
                push!(vx,vx4(t))
                push!(vy,vy4(t))
            elseif t4 < t
                push!(vx,vx5(t))
                push!(vy,vy5(t))
                if y5(t) < 0
                    break
                end
            end

            if floor(Int,cnt*dt*1000)>tr
                tr += 1
                if track
                    println("$(tr) ms")
                end
            end

            cnt += 1
        end

        return vx , vy , T , [vx1(t1) , vx2(t2) , vx3(t3) , vx4(t4)] , [vy1(t1) , vy2(t2) , vy3(t3) , vy4(t4)] , [t1 , t2 , t3 , t4]
    end
end

function regime2(dt , print::Bool , track::Bool , PorV::Bool)
    #To do Position input true and velocity input false for "PorV"

    temp_u2(t) = exp((c1/m+c2*v2_p/(6*m))*t)*((abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))^(2/3))
    du2(t) = v2_m/temp_u2(t)-(g/temp_u2(t))*numInt(temp_u2 , 0 , t)
    u2(t) = numInt(du2 , 0 , t)
    #du2(t) = lambda_p*exp(-phi_p*t)-g/phi_p
    dv2(t) = (v2_p/3)+(2*m*phi_p/(3*c2))*((2*kappa_p-eta_p*exp(-phi_p*t))*real(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))-imag(conj(r)*(1+2*im*(kappa_p+mu_p))*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p*exp(-phi_p*t))))/imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))
    #u2(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
    v2(t) = (v2_p/3)*t+(4*m/(3*c2))*log(abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))
    x2(t) = (v2(t)-u2(t))/2
    y2(t) = (v2(t)+u2(t))/2
    vx2(t) = (dv2(t)-du2(t))/2
    vy2(t) = (dv2(t)+du2(t))/2
    ryx2(t) = vy2(t)/vx2(t)-q

    temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*t)*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*t)))
    temp_vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
    vy3(t) = v3_y/temp_y3(t-t2)-(g/temp_y3(t-t2))*numInt(temp_y3 , 0 , t-t2)
    y3(t) = numInt(temp_vy3 , 0 , t-t2) + d3_y
    vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t-t2)))*real(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))
    #vy3(t) = gammay*exp(-omega_0*(t-t2))-g/omega_0
    x3(t) = (v3_x/2)*(t-t2)+(m/c2)*log(abs(imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))))+d3_x
    #y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t-t2)))-(g/omega_0)*(t-t2)+d3_y
    ryx3(t) = -vy3(t)/vx3(t)-q

    temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*((abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3))
    dv4(t) = v4_p/temp_v4(t)-(g/temp_v4(t))*numInt(temp_v4 , 0 , t)
    v4(t) = numInt(dv4 , 0 , t)
    #dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
    du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
    #v4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
    u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))
    x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
    y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
    vx4(t) = (dv4(t-t3)-du4(t-t3))/2
    vy4(t) = (dv4(t-t3)+du4(t-t3))/2
    rxy4(t) = -vx4(t)/vy4(t)-q

    vx5(t) = chi_m*exp(-c1*(t-t4)/(2*c2))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
    vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
    x5(t) = numInt(vx5 , t4 , t) + d5_x
    y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(abs(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega)))+d5_y

    q = .677269

    v2_p = v0*(sinpi(θ/2)+cospi(θ/2))
    v2_m = v0*(sinpi(θ/2)-cospi(θ/2))
    phi_p = (2*c1+c2*v2_p)/(2*m)
    lambda_p = v2_m + g/phi_p
    eta_p = (3*c2*lambda_p)/(2*m*phi_p)
    kappa_p = (3*c2*g)/(4*m*phi_p^2)
    mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
    r = -(1/(2*eta_p*mu_p))*(rWhittaker(im*kappa_p,im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p))
    
    t2 = forwardEval(ryx2 , 0 , .1)
    
    d3_x = x2(t2)
    d3_y = y2(t2)
    v3_x = vx2(t2)
    v3_y = vy2(t2)
    omega_0 = (c1+c2*v3_x)/m
    gammay = v3_y+g/omega_0
    xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
    delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
    epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
    p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))
    
    t3 = forwardEval(ryx3 , t2 , .1)
    
    d4_x = x3(t3)
    d4_y = y3(t3)
    v4_p = vy3(t3)+vx3(t3)
    v4_m = vy3(t3)-vx3(t3)
    phi_m = (2*c1-c2*v4_m)/(2*m)
    lambda_m = v4_p + g/phi_m
    eta_m = (3*c2*lambda_m)/(2*m*phi_m)
    kappa_m = (3*c2*g)/(4*m*phi_m^2)
    mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
    s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

    t4 = forwardEval(rxy4 , t3 , .1)
    
    d5_x = x4(t4)
    d5_y = y4(t4)
    v5_y = vy4(t4)
    chi_m = vx4(t4)
    omega_m = (c1-c2*v5_y)/m
    xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
    Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
    l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
    l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

    if print
        println("")

        println("c1 = ",c1)
        println("c2 = ",c2)
        println("t2 = ",t2)
        if y3(t3) > 0
            println("t3 = ",t3)
        end
        if y4(t4) > 0
            println("t4 = ",t4)
        end
        
        println("")

        println("phi_p = ",phi_p)
        println("lambda_p = ",lambda_p)
        println("eta_p = ",eta_p)
        println("kappa_p = ",kappa_p)
        println("mu_p = ",mu_p)
        println("r = ",r)

        println("")
        println("omega_0 = ", omega_0)
        println("gammay = ", gammay)
        println("xi_0 = ", xi_0)
        println("delta = ", delta)
        println("epsilon = ", epsilon)
        println("p = ",p)
        println("")

        if y3(t3) > 0
            println("phi_m = ", phi_m)
            println("lambda_m = ", lambda_m)
            println("eta_m = ", eta_m)
            println("kappa_m = ", kappa_m)
            println("mu_m = ", mu_m)
            println("s = ",s)
            println("")
        end
        
        if y4(t4) > 0
            println("chi_m = ", chi_m)
            println("omega_m = ", omega_m)
            println("xi_m = ", xi_m)
            println("Omega = ", Omega)
            println("l_p = ",l_p)
            println("l_m = ",l_m)
            println("")
        end
    end

    if PorV
        x = []
        y = []
        T = []
        tr = 0
        cnt = 0
        while true

            t = cnt*dt
            push!(T,t)
            
            if t <= t2
                push!(x,x2(t))
                push!(y,y2(t))
            elseif t2 < t <= t3
                push!(x,x3(t))
                push!(y,y3(t))
                if y3(t) < 0
                    break
                end
            elseif t3 < t <= t4
                push!(x,x4(t))
                push!(y,y4(t))
                if y4(t) < 0
                    break
                end
            elseif t4 < t
                push!(x,x5(t))
                push!(y,y5(t))
                if y5(t) < 0
                    break
                end
            end

            
            if floor(Int,cnt*dt*1000)>tr
                tr += 1
                if track
                    println("$(tr) ms")
                end
            end

            cnt += 1
        
        end

        xt = []
        yt = []
        Tt = []
        if y2(t2) > 0
            push!(xt , x2(t2))
            push!(yt , y2(t2))
            push!(Tt , t2)
        end
        if y3(t3) > 0
            push!(xt , x3(t3))
            push!(yt , y3(t3))
            push!(Tt , t3)
        end
        if y4(t4) > 0
            push!(xt , x4(t4))
            push!(yt , y4(t4))
            push!(Tt , t4)
        end

        return x , y , T , xt , yt , Tt
    else
        vx = []
        vy = []
        T = []
        tr = 0
        cnt = 0
        while true

            t = cnt*dt
            push!(T,t)

            if t <= t1
                push!(vx,vx1(t))
                push!(vy,vy1(t))
            elseif t1 < t <= t2
                push!(vx,vx2(t))
                push!(vy,vy2(t))
            elseif t2 < t <= t3
                push!(vx,vx3(t))
                push!(vy,vy3(t))
                if y3(t) < 0
                    break
                end
            elseif t3 < t <= t4
                push!(vx,vx4(t))
                push!(vy,vy4(t))
                if y4(t) < 0
                    break
                end
            elseif t4 < t
                push!(vx,vx5(t))
                push!(vy,vy5(t))
                if y5(t) < 0
                    break
                end
            end

            if floor(Int,cnt*dt*1000)>tr
                tr += 1
                if track
                    println("$(tr) ms")
                end
            end

            cnt += 1
        end

        xt = []
        yt = []
        Tt = []
        if y2(t2) > 0
            push!(xt , vx2(t2))
            push!(yt , vy2(t2))
            push!(Tt , t2)
        end

        if y3(t3) > 0
            push!(xt , vx3(t3))
            push!(yt , vy3(t3))
            push!(Tt , t3)
        end

        if y4(t4) > 0
            push!(xt , vx4(t4))
            push!(yt , vy4(t4))
            push!(Tt , t4)
        end

        return vx , vy , T , xt , yt , Tt
    end
end

function regime3(dt , print::Bool , track::Bool , PorV::Bool)
    #To do Position input true and velocity input false for "PorV"

    temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*t)*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*t)))
    temp_vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
    vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
    y3(t) = numInt(temp_vy3 , 0 , t)
    vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t)))*real(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t))))
    #vy3(t) = gammay*exp(-omega_0*(t))-g/omega_0
    x3(t) = (v3_x/2)*(t)+(m/c2)*log(abs(imag(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t))))))
    #y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t)))-(g/omega_0)*(t)
    ryx3(t) = -vy3(t)/vx3(t)-q

    temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*((abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3))
    dv4(t) = v4_p/temp_v4(t)-(g/temp_v4(t))*numInt(temp_v4 , 0 , t)
    v4(t) = numInt(dv4 , 0 , t)
    #dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
    du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
    #4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
    u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))
    x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
    y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
    vx4(t) = (dv4(t-t3)-du4(t-t3))/2
    vy4(t) = (dv4(t-t3)+du4(t-t3))/2
    rxy4(t) = -vx4(t)/vy4(t)-q

    vx5(t) = chi_m*exp(-c1*(t-t4)/(2*c2))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
    vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
    x5(t) = numInt(vx5 , t4 , t) + d5_x
    y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(abs(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega)))+d5_y

    q = .677269
    
    v3_x = v0*cospi(θ/2)
    v3_y = v0*sinpi(θ/2)
    omega_0 = (c1+c2*v3_x)/m
    gammay = v3_y+g/omega_0
    xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
    delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
    epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
    p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))
    
    t3 = forwardEval(ryx3 , 0 , .1)
    
    d4_x = x3(t3)
    d4_y = y3(t3)
    v4_p = vy3(t3)+vx3(t3)
    v4_m = vy3(t3)-vx3(t3)
    phi_m = (2*c1-c2*v4_m)/(2*m)
    lambda_m = v4_p + g/phi_m
    eta_m = (3*c2*lambda_m)/(2*m*phi_m)
    kappa_m = (3*c2*g)/(4*m*phi_m^2)
    mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
    s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

    t4 = forwardEval(rxy4 , t3 , .1)
    
    d5_x = x4(t4)
    d5_y = y4(t4)
    v5_y = vy4(t4)
    chi_m = vx4(t4)
    omega_m = (c1-c2*v5_y)/m
    xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
    Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
    l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
    l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

    if print
        println("")

        println("c1 = ",c1)
        println("c2 = ",c2)
        if y3(t3) > 0
            println("t3 = ",t3)
        end
        if y4(t4) > 0
        println("t4 = ",t4)
        end
        
        println("")
        println("omega_0 = ", omega_0)
        println("gammay = ", gammay)
        println("xi_0 = ", xi_0)
        println("delta = ", delta)
        println("epsilon = ", epsilon)
        println("p = ",p)
        println("")

        if y3(t3) > 0
            println("phi_m = ", phi_m)
            println("lambda_m = ", lambda_m)
            println("eta_m = ", eta_m)
            println("kappa_m = ", kappa_m)
            println("mu_m = ", mu_m)
            println("s = ",s)
            println("")
        end

        if y4(t4) > 0
            println("chi_m = ", chi_m)
            println("omega_m = ", omega_m)
            println("xi_m = ", xi_m)
            println("Omega = ", Omega)
            println("l_p = ",l_p)
            println("l_m = ",l_m)
            println("")
        end
    end

    if PorV
        x = []
        y = []
        T = []
        tr = 0
        cnt = 0
        while true

            t = cnt*dt
            push!(T,t)
            
            if t <= t3
                push!(x,x3(t))
                push!(y,y3(t))
                if y3(t) < 0
                    break
                end
            elseif t3 < t <= t4
                push!(x,x4(t))
                push!(y,y4(t))
                if y4(t) < 0
                    break
                end
            elseif t4 < t
                push!(x,x5(t))
                push!(y,y5(t))
                if y5(t) < 0
                    break
                end
            end

            
            if floor(Int,cnt*dt*1000)>tr
                tr += 1
                if track
                    println("$(tr) ms")
                end
            end

            cnt += 1
        
        end

        xt = []
        yt = []
        Tt = []
        
        if y3(t3) > 0
            push!(xt , x3(t3))
            push!(yt , y3(t3))
            push!(Tt , t3)
        end

        if y4(t4) > 0
            push!(xt , x4(t4))
            push!(yt , y4(t4))
            push!(Tt , t4)
        end

        return x , y , T , xt , yt , Tt
    else
        vx = []
        vy = []
        T = []
        tr = 0
        cnt = 0
        while true

            t = cnt*dt
            push!(T,t)

            if t <= t1
                push!(vx,vx1(t))
                push!(vy,vy1(t))
            elseif t1 < t <= t2
                push!(vx,vx2(t))
                push!(vy,vy2(t))
            elseif t2 < t <= t3
                push!(vx,vx3(t))
                push!(vy,vy3(t))
                if y3(t) < 0
                    break
                end
            elseif t3 < t <= t4
                push!(vx,vx4(t))
                push!(vy,vy4(t))
                if y4(t) < 0
                    break
                end
            elseif t4 < t
                push!(vx,vx5(t))
                push!(vy,vy5(t))
                if y5(t) < 0
                    break
                end
            end

            if floor(Int,cnt*dt*1000)>tr
                tr += 1
                if track
                    println("$(tr) ms")
                end
            end

            cnt += 1
        end

        xt = []
        yt = []
        Tt = []
        if y3(t3) > 0
            push!(xt , vx3(t3))
            push!(yt , vy3(t3))
            push!(Tt , t3)
        end

        if y4(t4) > 0
            push!(xt , vx4(t4))
            push!(yt , vy4(t4))
            push!(Tt , t4)
        end

        return vx , vy , T , xt , yt , Tt
    end
end