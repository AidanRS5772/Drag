using Plots
using SpecialFunctions
using QuadGK
using TickTock
using Printf

#=
For the functions below they calculate using a quadrature integral repersentations 
of function associated with the Bessel Function of the First Kind 
    x : input of the function (Real)
    o : digits of accuracy (Integer)
=#
function Bessel(x,z)
    #if b is False then the imaginary gamma values are taken if True the Rea
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

function quadDrag()
    T = floor(Int,rng/dt)-1

    y = [0.0]
    x = [0.0]

    vy=v0*sin(th)
    vx=v0*cos(th)

    for i = 1:T
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

        push!(x,nx)
        push!(y,ny)
    end

    return x,y
end

function y1(t)
    arg1 = alpha*real(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    arg2 = beta*imag(Bessel(xi*exp(-((c1+c2*am)/m)*t),im*zeta))
    return -(c1/(2*c2))*t + (m/c2)log(arg1-arg2)
end

function y2(t)
    return lam2-(lam1*m/c1)*exp(-c1*t/m)-g*m*t/c1
end

function y3(t)
    arg1 = mu*Bessel(psi*exp(-(c1-c2*ap)*t/m),omega)
    arg2 = nu*Bessel(psi*exp(-(c1-c2*ap)*t/m),-omega)
    return c1*t/(2*c2)-(m/c2)*log(abs(arg1-arg2))+log(abs(π*psi/(sin(π*omega))))+Lam
end

function x(t)
    if t < tau1
        return (m*v0*cos(th)/(c1+am*c2))*(1 - exp.(-(c1+am*c2)*t/m))
    elseif tau1 < t && t<tau2
        return del2-(del1*m/c1)*exp(-c1*t/m)
    else
        return -(sqrt(2)*m*psi/c2)*exp(-(c1-c2*ap)*t/m)+d
    end
end

function y(t)
     return y1(t)
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



c1 = linC*D
c2 = quadC*D^2

#=
Linear aproximations of y in the regime of 1 and 2 where "am" 
coresponds to regime 1 and "ap" to regime 3 
=#
am = v0*sin(th)
ap = -3

tp = sqrt(m/(g*c2))*atan(v0*sin(th)*sqrt(c2/(m*g)))
yp = (m/(2*c2))*log((v0^2)*((sin(th))^2)*c2/(m*g)+1)
tau1 = tp-.05
tau2 = tp+.05

println("\nPeak Time:")
println(tp)
println("\nPeak Height:")
println(yp)

global zeta = sqrt(4*g*m*c2-c1^2)/(2*(c1+c2*am))
global omega = sqrt(4*g*m*c2+c1^2)/(2*(c1+c2*am))
global xi = (c2*v0*cos(th))/(sqrt(2)*(c1+am*c2))
global pxi = xi*exp(-((c1+c2*am)/m)*tau1)
global psi = (c2*v0*cos(th)/(sqrt(2)*(c1-c2*ap)))*exp(-c2*(am*tau1+ap*tau2)/m)
global ppsi = (c2*v0*cos(th)/(sqrt(2)*(c1-c2*ap)))*exp(-(c1*tau2+c2*am*tau1)/m)

global var1 = sqrt(8)*tan(th)+((c1*sqrt(2))/(v0*c2))*sec(th)
global Xi  = real(Bessel(xi,im*zeta))*imag(Bessel(xi,-im*zeta+1)-Bessel(psi,-im*zeta-1))-imag(Bessel(psi,-im*zeta))*real(Bessel(psi,omega-1)-Bessel(psi,omega+1))
global alpha = (imag(Bessel(xi,im*zeta))*var1-imag(Bessel(xi,im*zeta+1) - Bessel(xi,im*zeta-1)))/Xi
global beta = (real(Bessel(xi,im*zeta))*var1-real(Bessel(xi,im*zeta+1) - Bessel(xi,im*zeta-1)))/Xi

global var2 = (alpha*real(Bessel(pxi,im*zeta+1) - Bessel(pxi,im*zeta-1))-beta*imag(Bessel(pxi,im*zeta+1) - Bessel(pxi,im*zeta-1)))/(alpha*real(Bessel(pxi,im*zeta))-beta*imag(Bessel(pxi,im*zeta)))
global del1 = v0*cos(th)*exp(-c2*am*tau1/m)
global del2 = m*v0*cos(th)*(1/(c1+c2*am)+exp(-(c1+c2*am)*tau1/m)*(1/c1-1/(c1+c2*am)))
global lam1 = (g*m/c1-c1/(2*c2))*exp(c1*tau1/m)+(v0*cos(th)/sqrt(8))*exp(-c2*am*tau1/m)*var2
global lam2 = (g*m/c1-c1/(2*c2))*(tau1+m/c1)+(m*v0*cos(th)/(sqrt(8)*c1))*exp(-(c1+c2*am)*tau1/m)*var2+(m/c2)*log(alpha*real(Bessel(pxi,im*zeta))-beta*imag(Bessel(pxi,im*zeta)))+(m/c2)*log((π*xi)/sinh(π*zeta))

global d = m*v0*cos(th)*(exp(-(c1+c2*am)*tau1/m)*(1/c1-1/(c1+c2*am))-exp(-(c1*tau2+c2*am*tau1)/m)*(1/c1-1/(c1-c2*ap))+1/(c1+c2*am))
global Lam = (lam1*c1/c2)*exp(-c1*tau2/m)+(g*c2/c1+c1/(2*m))*tau2-lam2
global mu = real(Bessel(ppsi,-omega)*(c1/(2*c2)+g*m/c1-lam1*exp(-c1*tau2/m))*(2*c2/(psi*(c1-c2*ap)))-Bessel(ppsi,-omega+1)+Bessel(ppsi,-omega-1))
global nu = real(Bessel(ppsi,omega)*(c1/(2*c2)+g*m/c1-lam1*exp(-c1*tau2/m))*(2*c2/(psi*(c1-c2*ap)))-Bessel(ppsi,omega+1)+Bessel(ppsi,omega-1))

println("\nzeta:")
println(zeta)

println("\nomega:")
println(omega)

println("\nxi:")
println(xi)

println("\npxi:")
println(pxi)

println("\npsi:")
println(psi)

println("\nppsi:")
println(ppsi)

println("\nXi:")
println(Xi)

println("\nalpha:")
println(alpha)

println("\nbeta:")
println(beta)

println("\ndel1:")
println(del1)

println("\ndel2:")
println(del2)

println("\nlam1:")
println(lam1)

println("\nlam2:")
println(lam2)

println("\nd:")
println(d)

println("\nLam:")
println(Lam)

println("\nmu:")
println(mu)

println("\nnu:")
println(nu)

#Length of time being utilised
global rng = tp*2
#Time step
global dt = .0001
#scale factor for Plot

n = floor(Int,rng/dt)
t = LinRange(0,rng,n)

vacY = -(g/2).*t.^2+sin(th)*v0.*t
vacX = cos(th)*v0.*t

linY = (m/c1)*(g*m/c1+v0*sin(th)).*(1 .- exp.(((-1)*c1/m).*t)) .- (g*m/c1).*t
linX = (m*v0*cos(th)/c1).*(1 .- exp.(((-1)*c1/m).*t))

quadX = quadDrag()[1]
quadY = quadDrag()[2]

aproxX = []
aproxY = []


track = 0
tick()
for i in 1:n
    push!(aproxX,x(t[i]))
    push!(aproxY,y(t[i]))
    if floor(Int,100*i/n) > track
        global track += 1
        @printf("%.f%% \n",100*i/n)
    end
end
tock()

#=
abserrorx = convert.(Float64,abs.(quadX-aproxX))
abserrory = convert.(Float64,abs.(quadY-aproxY))

relerrorx = convert.(Float64,abs.(quadX-aproxX)./abs.(aproxX))
relerrory = convert.(Float64,abs.(quadY-aproxY)./abs.(aproxY))

abserrorx = filter(!isnan,abserrorx)
abserrory = filter(!isnan,abserrory)
relerrorx = filter(!isnan,relerrorx)
relerrory = filter(!isnan,relerrory)

println("\n\nAbsolute Error X:")
println(maximum(abserrorx))
println("\n\nAbsolute Error Y:")
println(maximum(abserrory))
println("\n\nReletive Error X:")
println(maximum(relerrorx))
println("\n\nReletive Error Y:")
println(maximum(relerrory))
=#

p = plot(aproxX,aproxY)
#p = plot(quadX,quadY)
display(p)

