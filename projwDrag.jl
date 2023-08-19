using SpecialFunctions
using Cubature
using Printf

function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting, this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*a-(c2/m)*a*sqrt(a^2+b^2)
end

function qDSim(;dt = 2.0^-16,track = true)
    if track
        println("\nSimulation:\n")
    end

    time = 0
    t = []
    y = []
    x = []
    vY = []
    vX = []

    vy=v0*sinpi(θ/2)
    vx=v0*cospi(θ/2)

    ny = 0
    nx = 0

    while ny > -eps(Float64) 

        time += dt
        push!(t , time)

        push!(x,nx)
        push!(y,ny)
        push!(vX,vx)
        push!(vY,vy)


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

    return x , y , t , vX , vY
end

function secantRF(f::Function , init)
    tol = 1e-8
    
    p0 = init
    p1 = init + 2*tol
    p2 = 0

    fp0 = f(p0)
    fp1 = f(p1)

    while abs(p1-p0) > tol
        p2 = p1 - fp1*((p1-p0)/(fp1-fp0))
        
        p0 = p1
        p1 = p2

        fp0 = fp1
        fp1 = f(p2)
    end

    return p2
end

function d1Int(f::Function , t)
    out , _ = hquadrature(f , 0 , t)
    return out
end

function d2Int(f::Function , (x1 , x2))
    out , _ = hcubature(f , (0,0), (x1,x2))
    return out
end

function triTorect(f1)
    return function f2(x)
        return x[1] * f1(x[1] * x[2]) / f1(x[1])
    end
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
    n = 0

    while true
        val = (((-1)^n)/(GammaAll(n+z+1)*GammaAll(n+1)))*(x/2)^(2*n)
        sum += val
        if (abs(val)/abs(sum) < 1e-15) break end
        n += 1
    end
    if isa(x , Complex) || isa(z , Complex)
        return convert(ComplexF64 , ((x/2)^z)*sum)
    else
        return ((x/2)^z)*sum
    end
end

function rWhittaker(a,b,x)
    sum = 0
    n = 0

    while true
        val = (GammaAll(b-a+1/2+n)/GammaAll(2*b+1+n))*(big(x)^n/GammaAll(n+1))
        sum += val
        if (abs(val)/abs(sum) < 1e-10) break end
        n += 1
    end

    return convert(ComplexF64 , sqrt(x)*exp(-x/2)*(GammaAll(2*b+1)/GammaAll(b-a+1/2))*sum)
end

#Regime 1

#old x1
x1(t) = (chi_p / omega_p) * (1 - exp(-omega_p * t))

#new x1
#temp_x1(t) = 

y1(t) = -(c1/(2*c2))*t+(m/c2)*log(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))
vx1(t) = chi_p*exp(-omega_p*t)
vy1(t) = -(c1/(2*c2))+(m*xi_p*omega_p/(2*c2))*exp(-omega_p*t)*((imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta+1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))-(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta-1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))
rxy1(t) = vx1(t)-vy1(t)*q

u2(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
v2(t) = (v2_p/3)*t+(4*m/(3*c2))*log(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t))))
du2(t) = lambda_p*exp(-phi_p*t) - g/phi_p
dv2(t) = (v2_p/3)+(2*m*phi_p/(3*c2))*((2*kappa_p-eta_p*exp(-phi_p*t))*real(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))-imag(conj(r)*(1+2*im*(kappa_p+mu_p))*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p*exp(-phi_p*t))))/imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))
x2(t) = (v2(t-t1)-u2(t-t1))/2 + d2_x
y2(t) = (v2(t-t1)+u2(t-t1))/2 + d2_y
vx2(t) = (dv2(t-t1)-du2(t-t1))/2
vy2(t) = (dv2(t-t1)+du2(t-t1))/2
ryx2(t) = vy2(t)-vx2(t)*q

x3(t) = (v3_x/2)*(t-t2)+(m/c2)*log(imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))+d3_x
y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t)))-(g/omega_0)*(t)
vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t-t2)))*real(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))
vy3(t) = gammay*exp(-omega_0*(t))-g/omega_0
ryx3(t) = vy3(t)+vx3(t)*q

v4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t))))
dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
vx4(t) = (dv4(t-t3)-du4(t-t3))/2
vy4(t) = (dv4(t-t3)+du4(t-t3))/2
rxy4(t) = vx4(t)+vy4(t)*q

x5(t) = (chi_m/omega_m)*(1-exp(-omega_m*(t)))
y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(l_p*Bessel(xi_m*exp(-omega_m*(t)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t)),-Omega))+d5_y
vx5(t) = chi_m*exp(-omega_m*(t))
vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))

function qDAproxP(t)
    if (t <= t1) && (t1 != 0)
        return x1(t) , y1(t)
    elseif (t <= t2) && (t2 != 0)
        return x2(t) , y2(t)
    elseif t <= t3
        return x3(t) , y3(t)
    elseif t <= t4
        return x4(t) , y4(t)
    else
        return x5(t) , y5(t)
    end
end

function qDAproxV(t)
    if (t <= t1) && (t1 != 0)
        return vx1(t) , vy1(t)
    elseif (t <= t2) && (t2 != 0)
        return vx2(t) , vy2(t)
    elseif t <= t3
        return vx3(t) , vy3(t)
    elseif t <= t4
        return vx4(t) , vy4(t)
    else
        return vx5(t) , vy5(t)
    end
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
    global q = .677269
end

function d1Vals(vx , vy , print)
    global chi_p = vx
    global omega_p = (c1 + c2 * vy) / m
    global xi_p = (c2 * chi_p) / (sqrt(2) * omega_p * m)
    global zeta = sqrt(4*g*c2*m - c1^2) / (2*m*omega_p)
    global k = -(π/(omega_p*sinh(π*zeta)))*(Bessel(xi_p , im*zeta)*((c2*v0*sinpi(θ/2)/m)+(c1/(2*m))-im*zeta*omega_p)+xi_p*omega_p*Bessel(xi_p,im*zeta-1))

    if print
        println("chi_p = ",chi_p)
        println("omega_p = ",omega_p)
        println("xi_p = ",xi_p)
        println("zeta = ",zeta)
        println("k = ",k)
    end
end

function d2Vals(x , y , vx , vy , print)
    global d2_x = x
    global d2_y = y
    global v2_p = vy+vx
    global v2_m = vy-vx
    global phi_p = (2*c1+c2*v2_p)/(2*m)
    global lambda_p = v2_m + g/phi_p
    global eta_p = (3*c2*lambda_p)/(2*m*phi_p)
    global kappa_p = (3*c2*g)/(4*m*phi_p^2)
    global mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
    global r = -(1/(2*eta_p*mu_p))*(rWhittaker(im*kappa_p,im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p))

    if print
        println("phi_p = ",phi_p)
        println("lambda_p = ",lambda_p)
        println("eta_p = ",eta_p)
        println("kappa_p = ",kappa_p)
        println("mu_p = ",mu_p)
        println("r = ",r)
    end
end

function d3Vals(x , y , vx , vy , print)
    global d3_x = x
    global d3_y = y
    global v3_x = vx
    global v3_y = vy
    global omega_0 = (c1+c2*v3_x)/m
    global gammay = v3_y+g/omega_0
    global xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
    global delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
    global epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
    global p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))

    if print
        println("omega_0 = ", omega_0)
        println("gammay = ", gammay)
        println("xi_0 = ", xi_0)
        println("delta = ", delta)
        println("epsilon = ", epsilon)
        println("p = ", p)
    end
end

    
function d4Vals(x , y , vx , vy , print)
    global d4_x = x
    global d4_y = y
    global v4_p = vy+vx
    global v4_m = vy-vx
    global phi_m = (2*c1-c2*v4_m)/(2*m)
    global lambda_m = v4_p + g/phi_m
    global eta_m = (3*c2*lambda_m)/(2*m*phi_m)
    global kappa_m = (3*c2*g)/(4*m*phi_m^2)
    global mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
    global s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))
    
    if print
        println("phi_m = ", phi_m)
        println("lambda_m = ", lambda_m)
        println("eta_m = ", eta_m)
        println("kappa_m = ", kappa_m)
        println("mu_m = ", mu_m)
        println("s = ",s)
    end
end

function d5Vals(x , y , vx , vy , print)
    global d5_x = x
    global d5_y = y
    global v5_y = vy
    global chi_m = vx
    global omega_m = (c1-c2*v5_y)/m
    global xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
    global Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
    global l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
    global l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

    if print
        println("chi_m = ", chi_m)
        println("omega_m = ", omega_m)
        println("xi_m = ", xi_m)
        println("Omega = ", Omega)
        println("l_p = ",l_p)
        println("l_m = ",l_m)
    end
end


function preCalc(;print = true)
    push = .01

    if 1 > θ >= acot(q)*2/π
        d1Vals(v0*cospi(θ/2) , v0*sinpi(θ/2) , print)
        global t1 = secantRF(rxy1 , 0)

        if print
            println("t1 = ",t1)
            println("")
        end

        d2Vals(x1(t1) , y1(t1) , vx1(t1) , vy1(t1) , print)
        global t2 = secantRF(ryx2 , t1)

        if print
            println("t2 = ",t2)
            println("")
        end

        d3Vals(x2(t2) , y2(t2) , vx2(t2) , vy2(t2) , print)
        global t3 = secantRF(ryx3 , t2)

        y3_t3 = y3(t3)
        if y3_t3 >= 0
            if print
                println("t3 = ",t3)
                println("")
            end

            d4Vals(x3(t3) , y3_t3 , vx3(t3) , vy3(t3) , print)
            global t4 = secantRF(rxy4 , t3)
        end

        y4_t4 = y4(t4)
        if y4_t4 >= 0
            if print
                println("t4 = ",t4)
                println("")
            end

            d5Vals(x4(t4) , y4_t4 , vx4(t4) , vy4(t4) , print)
        end

    elseif acot(q)*2/π >= θ > atan(q)*2/π

        d2Vals(0 , 0 , v0*cospi(θ/2) , v0*sinpi(θ/2) , print)
        global t2 = secantRF(ryx2 , 0)

        if print
            println("t2 = ",t2)
            println("")
        end
        
        d3Vals(x2(t2) , y2(t2) , vx2(t2) , vy2(t2) , print)
        global t3 = secantRF(ryx3 , t2)

        y3_t3 = y3(t3)
        if y3_t3 >= 0
            if print
                println("t3 = ",t3)
                println("")
            end

            d4Vals(x3(t3) , y3_t3 , vx3(t3) , vy3(t3) , print)
            global t4 = secantRF(rxy4 , t3)
        end

        y4_t4 = y4(t4)
        if y4_t4 >= 0
            if print
                println("t4 = ",t4)
                println("")
            end

            d5Vals(x4(t4) , y4_t4 , vx4(t4) , vy4(t4) , print)
        end

        global t1 = 0

    elseif atan(q)*2/π >= θ > 0

        d3Vals(x2(t2) , y2(t2) , vx2(t2) , vy2(t2) , print)
        global t3 = secantRF(ryx3 , 0)

        y3_t3 = y3(t3)
        if y3_t3 >= 0
            if print
                println("t3 = ",t3)
                println("")
            end

            d4Vals(x3(t3) , y3_t3 , vx3(t3) , vy3(t3) , print)
            global t4 = secantRF(rxy4 , t3)
        end

        y4_t4 = y4(t4)
        if y4_t4 >= 0
            if print
                println("t4 = ",t4)
                println("")
            end

            d5Vals(x4(t4) , y4_t4 , vx4(t4) , vy4(t4) , print)
        end

        global t1 = 0
        global t2 = 0
    end
end
