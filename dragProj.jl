using Plots
using SpecialFunctions
using QuadGK
using TickTock
using Printf


function Bessel(x,z)
    n=100
    sum = 0
    term  = (x/2)^z
    for i in 0:n
        term *= (-1)*x^2/(n*4)
        sum += term/gamma(n+z+1)
    end
    return convert(ComplexF64,sum)
end

function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*a-(c2/m)*a*sqrt(a^2+b^2)
end

function quadDragSim()
    y = []
    x = []

    vy=v0*sin(th)
    vx=v0*cos(th)

    ny = 0
    nx = 0

    while ny > -.001

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

    return x,y
end

function y1()

end

##################################
g = 9.8
th = (pi/2)*.9
v0 = 3
m = 1

linC = .00016
quadC = .25
D =.05
#################################

global c1 = linC*D
global c2 = quadC*D^2

simX = quadDragSim()[1]
simY = quadDragSim()[2]

global dt = .0001



