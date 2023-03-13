using Plots
using SpecialFunctions
using QuadGK
using TickTock
using Printf

function alt(n)
    if n%2 == 0
        return 1
    else 
        return -1
    end
end

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

function quadDragSim()
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
    return -(c1*t/(2*c2))+(m/c2)*log(arg1-arg2)-(m/c2)*log(Xi)
end

function x1(t)
    v1 = m*v0*cos(th)/(c1+c2*am)
    v2 = (c1+c2*am)/m
    return v1*(1-exp(-v2*t))
end

function quadDragAprox(T)
    x = []
    y = []

    for t in T
        push!(x,x1(t))
        push!(y,y1(t))
    end

    return x , y
end

##################################
g = 9.8
th = (pi/2)*.9
v0 = 2
m = 1

linC = .00016
quadC = .25
D =.05
#################################

global c1 = linC*D
global c2 = quadC*D^2
global am = v0*sin(th)
global zeta = sqrt(c2*g*m-(c1^2)/4)/(c1+c2*am)

global xi = c2*v0*cos(th)/(sqrt(2)*(c1+am*c2))
global Xi = imag(Bessel(xi,im*zeta))*real(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1))-real(Bessel(xi,im*zeta))*imag(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1))
global alpha = imag(Bessel(xi,im*zeta))*(2*sqrt(2)*tan(th)+(c1*sqrt(2)/(v0*c2))*sec(th))-imag(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1))
global beta = real(Bessel(xi,im*zeta))*(2*sqrt(2)*tan(th)+(c1*sqrt(2)/(v0*c2))*sec(th))-real(Bessel(xi,im*zeta+1)-Bessel(xi,im*zeta-1))

global dt = .001

simX = quadDragSim()[1]
simY = quadDragSim()[2]
n = quadDragSim()[3]

t = LinRange(0,dt*n,n)
aproxX = quadDragAprox(t)[1]
aproxY = quadDragAprox(t)[2]

#=
p1 = plot(simX,simY)
p2 = plot(aproxX,aproxY)
display(plot(p1,p2))
=#

abserrorX = abs.(simX .- aproxX)
abserrorY = abs.(simY .- aproxY)
relerrorX = abserrorX./simX
relerrorY = abserrorY./simY

p1 = plot(t,[abserrorX abserrorY],label = ["Absolute Error X" "Absolute Error Y"])
p2 = plot(t,[relerrorX relerrorY],label = ["Reletive Error X" "Reletive Error Y"])
display(plot(p1,p2))
