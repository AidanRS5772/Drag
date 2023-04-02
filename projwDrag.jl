using SpecialFunctions
using Printf

#=
function Bessel(x,z)
    sum = 0
    n=100
    for i in 0:n
        sum+= (((-1)^i)/(factorial(big(i))*gamma(z+i+1)))*(x/2)^(2*i+z)
    end
    return sum
end
=#

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

f(x,k) = x^3-(x^2)*(2*k+σ2)+x*(k^2+k*(σ2-1)-2)+(σ2+k)*(k+1)-1

function product(l)
    prod = 1
    for r in 1:l+1
        prod *= 1/(r^2+ζ^2)
    end
    return prod
end
 
function fFac(z,k)
    prod = 1
    for r in 0:k+1
        prod *= z+1-r
    end
    return prod
end

function y1Sum(t)
    sum1 = 0 
    n = 20
    for k in 0:n
        sum2 = 0
        for l in 0:k
            prod1 = product(l)
            prod2 = product(k-1)
            diff = cos(σ3*t)*imag(f(l+im*ζ,k)*fFac(l-im*ζ,k))-sin(σ3*t)*real(f(l-im*ζ,k)*fFac(l+im*ζ,k))
            sum2 += ((-1)^l)*(exp(-2*l*σ1*t)/(factorial(k-l)*factorial(l)))*prod1*prod2*diff
        end
        sum1 += ((ξ/4)^(k+1))*sum2
    end
    return sum1
end

y1(t) = -(c1/(2*c2))*t+(m/c2)*log(y1Sum(t))+(m/c2)*log(2/ζ^2)

function x1(t)
    return (sqrt(2)*m*ξ/c2)*(1-exp(-(c1+c2*am)*t/m))
end

function quadDragAprox(T)
    x = []
    y = []
    track = 0
    cnt = 0
    n = length(T)
    for t in T
        if t < τ
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
    global θ = theta
    global v0 = velocity
    global m = mass
    global D = diameter

    global g = 9.8
    global linC = .00016
    global quadC = .25
    global c1 = linC*D
    global c2 = quadC*D^2

    global ep = .1

    tp = quadDragSim(.0001)[5]
    global τ = tp + ep

    #free paremter
    global am = v0*sin(th)

    global ζ = sqrt(c2*g*m-(c1^2)/4)/(c1+c2*am)
    global ξ = c2*v0*cos(th)/(sqrt(2)*(c1+am*c2))
    global σ1 = (c1+c2*v0)/m
    global σ2 = 2*sqrt(2)*tan(θ)+c2*sqrt(2)*sec(θ)/(v0*c2)
    global σ3 = sqrt(c2*g/m-(c1/(2*m))^2)
    

    #=
    global ap = (y1(τ+eps(Float64))-y1(τ-eps(Float64)))/(2*eps(Float64))
    
    global Ω  = sqrt(c2*g*m+(c1^2)/4)/(c1-c2*ap)
    global ψ = (c2*v0*cos(th)/(sqrt(2)*(c1-c2*ap)))*exp(-(ap+am)*c2*τ/m)
    global pψ = (c2*v0*cos(th)/(sqrt(2)*(c1-c2*ap)))*exp(-(c1+c2*am)*τ/m)
    global pξ = (c2*v0*cos(th)/(sqrt(2)*(c1+c2*am)))*exp(-(c1+c2*am)*τ/m)
    global Δ = (α*real(Bessel(pξ,im*ζ+1)-Bessel(pξ,im*ζ-1))-β*imag(Bessel(pξ,im*ζ+1)-Bessel(pξ,im*ζ-1)))/(α*real(Bessel(pξ,im*ζ))-β*imag(Bessel(pξ,im*ζ)))
    global μ = Bessel(pψ,-Ω)*((2*sqrt(2)*c1/(c2*v0))*sec(th)*exp((c1+c2*am)*τ/m) + Δ)+Bessel(pψ,-Ω-1)-Bessel(pψ,-Ω+1)
    global ν = Bessel(pψ,Ω)*((2*sqrt(2)*c1/(c2*v0))*sec(th)*exp((c1+c2*am)*τ/m) + Δ)+Bessel(pψ,Ω-1)-Bessel(pψ,Ω+1)
    =#
end
