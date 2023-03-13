using Plots
using SpecialFunctions
using QuadGK
using TickTock
using Printf

function Bessel(x,z)
    sum = 0
    n=100
    for i in 0:n
        sum+= (((-1)^i)/(factorial(big(i))*gamma(i+z+1)))*(x/2)^(2*i+z)
    end
    return convert(ComplexF64,sum)
end


function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*a-(c2/m)*a*sqrt(a^2+b^2)
end

function quadDragSim(dt)
    cnt = 0

    y = []
    x = []

    vy=v0*sin(th)
    vx=v0*cos(th)

    ny = 0
    nx = 0

    while ny > -eps(Float64)

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

    return deleteat!(x,length(x)) , deleteat!(y,length(y)) , cnt-1
end

function y1(t)
    arg1 = alpha*real(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    arg2 = beta*imag(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    return -(c1*t/(2*c2))+(m/c2)*log(arg1-arg2)
end

function x1(t)
    v1 = m*v0*cos(th)/(c1+c2*am)
    v2 = (c1+c2*am)/m
    return v1*(1-exp(-v2*t))
end

function y2(t)
    return -lam1*(m/c1)*exp(-c1*t/m)-g*m*t/c1+lam2
end

function x2(t)
    return -del1*(m/c1)*exp(-c1*t/m)+del2
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
        else
            push!(x,x2(t))
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

function plotError(dt,cutOff)
    sim = quadDragSim(dt)

    simX = sim[1]
    simY = sim[2]
    n = sim[3]

    t = LinRange(0,dt*n,n)
    aprox = quadDragAprox(t)
    aproxX = aprox[1]
    aproxY = aprox[2]

    abserrorX = abs.(simX .- aproxX)
    abserrorY = abs.(simY .- aproxY)
    relerrorX = deleteat!(abserrorX./simX,floor(Int,n*cutOff):n)
    relerrorY = deleteat!(abserrorY./simY,floor(Int,n*cutOff):n)
    relt = deleteat!(collect(t),floor(Int,n*cutOff):n)

    ep1 = plot(t,[abserrorX abserrorY],legend = false)
    ep2 = plot(relt,[relerrorX relerrorY],legend = false)
    return ep1 , ep2
end

##################################
global g = 9.8
global th = (pi/2)*.9
global v0 = 2
global m = 1

global linC = .00016
global quadC = .25
global D =.05
#################################

global c1 = linC*D
global c2 = quadC*D^2
tp = sqrt(m/(g*c2))*atan(v0*sin(th)*sqrt(c2/(m*g)))
global tau1 = tp-.05
global tau2 = tp+.05
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

dt = .0001
sim = quadDragSim(dt)
simX = sim[1]
simY = sim[2]
n = sim[3]

t = LinRange(0,dt*n,n)
aprox = quadDragAprox(t)
aproxX = aprox[1]
aproxY = aprox[2]

p = plot(simX,simY,label = "Simulation")
display(plot(p,aproxX,aproxY,label = "Aproximation"))

#=
N = 3
absErrors = Array{Plots.Plot, 1}(undef, N);
relErrors = Array{Plots.Plot, 1}(undef, N);

for i in 1:N
    dt = 10.0^(-i-1)
    @printf("\n\n\n %.f \n\n\n",dt)
    absErrors[i] = plotError(dt,.8)[1]
    relErrors[i] = plotError(dt,.8)[2]
end
error = hcat(absErrors,relErrors)
display(plot(error...,layout = (2,N)))
=#