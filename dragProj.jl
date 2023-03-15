using Plots
using Interpolations
using SpecialFunctions
using TickTock
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

function quadDragSim()
    dt = 10.0^-6
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

    return x , y , time , cnt
end

function y1(t)
    arg1 = alpha*real(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    arg2 = beta*imag(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    return -(c1*big(t)/(2*c2))+(m/c2)*log(arg1-arg2)
end

function x1(t)
    v1 = m*v0*cos(th)/(c1+c2*am)
    v2 = (c1+c2*am)/m
    return v1*(1-exp(-v2*t))
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

function instVals()
    global c1 = linC*D
    global c2 = quadC*D^2
    tp = sqrt(m/(g*c2))*atan(v0*sin(th)*sqrt(c2/(m*g)))
    global tau1 = tp-ep1
    global tau2 = tp+ep2

    global am = v0*sin(th)
    global zeta = sqrt(c2*g*m-(c1^2)/4)/(c1+c2*am)
    global xi = c2*v0*cos(th)/(sqrt(2)*(c1+am*c2))
    global Xi = imag(Bessel(xi,im*zeta))*real(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1))-real(Bessel(xi,im*zeta))*imag(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1))
    global alpha = (imag(Bessel(xi,im*zeta))*(2*sqrt(2)*tan(th)+(c1*sqrt(2)/(v0*c2))*sec(th))-imag(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1)))/Xi
    global beta = (real(Bessel(xi,im*zeta))*(2*sqrt(2)*tan(th)+(c1*sqrt(2)/(v0*c2))*sec(th))-real(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1)))/Xi

    global pxi = xi*exp(-(c1+c2*am)*tau1/m)
    global var1 = (alpha*real(Bessel(pxi,im*zeta+1)-Bessel(pxi,im*zeta-1))-beta*imag(Bessel(pxi,im*zeta+1)-Bessel(pxi,im*zeta-1)))/(alpha*real(Bessel(pxi,im*zeta))-beta*imag(Bessel(pxi,im*zeta)))
    global del1 = v0*cos(th)*exp(-c2*am*tau1/m)
    global del2 = m*v0*cos(th)*(1/(c1+c2*am)+exp(-(c1+c2*am)*tau1/m)*(1/c1-1/(c1+c2*am)))
    global lam1 = (g*m/c1-c1/(2*c2))*exp(c1*tau1/m)+(v0*cos(th)/(2*sqrt(2)))*exp(-c2*am*tau1/m)*var1
    global lam2 = (g*m/c1-c1/(2*c2))*(tau1+m/c1)+(m*v0*cos(th)/(2*sqrt(2)*c1))*exp(-(c1+c2*am)*tau1/m)*var1+(m/c2)*log(alpha*real(Bessel(pxi,im*zeta))-beta*imag(Bessel(pxi,im*zeta)))

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

function plotError(dt,cutOff,simX,simY,time,cnt)

    t = LinRange(0,time,floor(Int,time/dt))
    aprox = quadDragAprox(t)
    aproxX = convert.(Float64 , aprox[1])
    aproxY = convert.(Float64 , aprox[2])

    splineX = CubicSplineInterpolation(t, aproxX)
    splineY = CubicSplineInterpolation(t, aproxY)
    splineT = big.(collect(LinRange(0,time,cnt)))
    linAproxX = big.(collect(splineX(splineT)))
    linAproxY = big.(collect(splineY(splineT)))

    abserrorX = abs.(simX .- linAproxX)
    abserrorY = abs.(simY .- linAproxY)
    relerrorX = abserrorX./simX
    relerrorY = abserrorY./simY
    relT = collect(splineT)

    nrelerrorX = relerrorX[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    nrelerrorY = relerrorY[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]
    nrelT = relT[floor(Int,cnt*cutOff):floor(Int,cnt*(1-cutOff))]

    ep1 = plot(splineT,[abserrorX abserrorY],legend = false)
    ep2 = plot(nrelT,[nrelerrorX nrelerrorY],legend = false)
    return ep1 , ep2
end

function manyErrorPlots(N,b,cutOff)
    absErrors = Array{Plots.Plot, 1}(undef, N);
    relErrors = Array{Plots.Plot, 1}(undef, N);

    instVals()

    sim = quadDragSim()

    simX = sim[1]
    simY = sim[2]
    time = sim[3]
    cnt = sim[4]

    for i in 1:N
        dt = 10.0^(-i-b+1)
        @printf("\n\n\n!!New Plot!!   dt:%.6f \n\n\n",dt)
        error = plotError(dt,cutOff,simX,simY,time,cnt)
        absErrors[i] = error[1]
        relErrors[i] = error[2]
    end
    error = hcat(absErrors,relErrors)
    return plot(error...,layout = (2,N))
end

function plotProj(dt)
    instVals()

    sim = quadDragSim()

    simX = sim[1]
    simY = sim[2]
    time = sim[3]

    t = LinRange(0,time,floor(Int,time/dt))
    aprox = quadDragAprox(t)
    aproxX = aprox[1]
    aproxY = aprox[2]

    p = plot(simX,simY,label = "Simulation")
    return plot(p,aproxX,aproxY,label = "Aproximation")
end

##################################
global g = 9.8
global th = (pi/2)*.9
global v0 = 2
global m = 1

global linC = .00016
global quadC = .25
global D =.05
global ep1 = .005
#certian ep2 values throw errors
global ep2 = .02
#################################