using SpecialFunctions
using Printf

function Bessel(x,z)
    sum = 0
    n=100
    for i in 0:n
        sum+= (((-1)^i)/(factorial(big(i))*gamma(i+z+1)))*(x/2)^(2*i+z)
    end
    return sum
end


function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*big(a)-(c2/m)*big(a)*sqrt(big(a)^2+big(b)^2)
end

function quadDragSim(dt = 10.0^-6)
    time = 0
    cnt = 0

    y = []
    x = []

    vy=v0*sin(th)
    vx=v0*cos(th)

    ny = 0
    nx = 0

    while ny > -eps(Float64) 

        time += dt
        cnt += 1

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

    return x , y , time , cnt , ymaxi*time/cnt
end

function y1(t)
    arg1 = alpha*real(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    arg2 = beta*imag(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    return -(c1*big(t)/(2*c2))+(m/c2)*log(arg1-arg2)+(m/c2)*log(π*xi/(2*sinh(π*zeta)))
end

function x1(t)
    return (sqrt(2)*m*xi/c2)*(1-exp(-(c1+c2*am)*t/m))
end

function y2(t)
    return -lam1*(m/c1)*exp(-c1*big(t)/m)-g*m*big(t)/c1+lam2
end

function x2(t)
    return -del1*(m/c1)*exp(-c1*t/m)+del2
end

function y3(t)
    arg1 = mu*Bessel(psi*exp(-(c1-c2*ap)*t/m),omega)
    arg2 = nu*Bessel(psi*exp(-(c1-c2*ap)*t/m),-omega)
    return c1*big(t)/(2*c2)-(m/c2)*log(arg1-arg2)+(m/c2)*(log(Psi)-Lam)
end

function x3(t)
    return -m*chi3/(c1-c2*ap)*exp(-(c1-c2*ap)*t/m)+d
end

function quadDragAprox(T)
    x = []
    y = []
    track = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < tau1
            push!(x,x1(t))
            push!(y,y1(t))
        elseif (t > tau1) && (t < tau2)
            push!(x,x1(t))
            push!(y,y1(t))
        else
            push!(x,x1(t))
            push!(y,y1(t))
        end

        cnt += 1
        if floor(Int,100*cnt/n) > track
            track += 1
            @printf("%.f%% \n",100*cnt/n)
        end
    end

    return x , y
end

function instInputs(;theta = (pi/2)*.95 , velocity = 1.0 , mass = 1.0 , diameter = .05)
    #input values
    global th = theta
    global v0 = velocity
    global m = mass
    global D = diameter

    global g = 9.8
    global linC = .00016
    global quadC = .25
    global c1 = linC*D
    global c2 = quadC*D^2

    #free parameters
    global ep1 = .05
    global ep2 = .05

    tp = quadDragSim(.0001)[5]
    global tau1 = tp - ep1
    global tau2 = tp + ep2

    #free paremter
    global am = v0*sin(th)

    global zeta = sqrt(c2*g*m-(c1^2)/4)/(c1+c2*am)
    global xi = c2*v0*cos(th)/(sqrt(2)*(c1+am*c2))
    global alpha = (imag(Bessel(xi,im*zeta))*(2*sqrt(2)*tan(th)+(c1*sqrt(2)/(v0*c2))*sec(th))-imag(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1)))
    global beta = (real(Bessel(xi,im*zeta))*(2*sqrt(2)*tan(th)+(c1*sqrt(2)/(v0*c2))*sec(th))-real(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1)))

    global pxi = xi*exp(-(c1+c2*am)*tau1/m)
    global var1 = (alpha*real(Bessel(pxi,im*zeta+1)-Bessel(pxi,im*zeta-1))-beta*imag(Bessel(pxi,im*zeta+1)-Bessel(pxi,im*zeta-1)))/(alpha*real(Bessel(pxi,im*zeta))-beta*imag(Bessel(pxi,im*zeta)))
    global del1 = v0*cos(th)*exp(-c2*am*tau1/m)
    global del2 = m*v0*cos(th)*(1/(c1+c2*am)+exp(-(c1+c2*am)*tau1/m)*(1/c1-1/(c1+c2*am)))
    global lam1 = (g*m/c1-c1/(2*c2))*exp(c1*tau1/m)+(v0*cos(th)/(2*sqrt(2)))*exp(-c2*am*tau1/m)*var1
    global lam2 = (g*m/c1-c1/(2*c2))*(tau1+m/c1)+(m*v0*cos(th)/(2*sqrt(2)*c1))*exp(-(c1+c2*am)*tau1/m)*var1+(m/c2)*log(alpha*real(Bessel(pxi,im*zeta))-beta*imag(Bessel(pxi,im*zeta)))

    #free parameter
    global ap = lam1*exp(-c1*tau2/m)-g*m/c1

    global omega = sqrt(c2*g*m+c1^2/4)/(c1-c2*ap)
    global chi3 = del1*exp(-c2*ap*tau2/m)
    global d = m*del1*exp(-c1*tau2/m)*(1/(c1-c2*ap)-1/c1)+del2
    global psi = (c2*v0*cos(th)/(sqrt(2)*(c1-c2*ap)))*exp(-(c2/m)*(am*tau1+ap*tau2))
    global ppsi = (c2*v0*cos(th)/(sqrt(2)*(c1-c2*ap)))*exp(-(c1*tau2+c2*am*tau1)/m)
    global Psi = Bessel(ppsi,-omega)*(Bessel(ppsi,omega+1)-Bessel(ppsi,omega-1))-Bessel(ppsi,omega)*(Bessel(ppsi,-omega+1)-Bessel(ppsi,-omega-1))
    global Lam = (lam1*c2/c1)*exp(-c1*tau2/m)+(g*c2/c1+c1/(2*m))tau2-c2*lam2/m
    global mu = Bessel(ppsi,-omega)*(c1/(2*c2)+g*m/c1-lam1*exp(-c1*tau2/m))*(2*c2/((c1-c2*ap)*psi))-Bessel(ppsi,-omega+1)+Bessel(ppsi,-omega-1)
    global nu = Bessel(ppsi,omega)*(c1/(2*c2)+g*m/c1-lam1*exp(-c1*tau2/m))*(2*c2/((c1-c2*ap)*psi))-Bessel(ppsi,omega+1)+Bessel(ppsi,omega-1)
end
