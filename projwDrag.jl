using SpecialFunctions
using Printf

function Bessel(x,z)
    sum = 0
    n=100
    for i in 0:n
        sum+= (((-1)^i)/(factorial(big(i))*gamma(z+i+1)))*(x/2)^(2*i+z)
    end
    return sum
end

function dragEq(a,b)
    #put as the first entry the variable you want the equation to be repersenting this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1/m)*big(a)-(c2/m)*big(a)*sqrt(big(a)^2+big(b)^2)
end

function quadDragSim(;dt = 10.0^-4)
    time = 0
    cnt = 0

    y = []
    x = []

    vy=v0*sin(θ)
    vx=v0*cos(θ)

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


function y1Sum(t)
    am = v0*sin(θ)
    ζ = sqrt(c2*g*m-(c1^2)/4)/(c1+c2*am)
    ξ = c2*v0*cos(θ)/(sqrt(2)*(c1+am*c2))
    ψ = 2*sqrt(2)*tan(θ)+c2*sqrt(2)*sec(θ)/(v0*c2)
    χ = (c1+c2*am)/m
    sum1 = 0 
    n = 20
    for k in 0:n
        sum2 = 0
        for l in 0:k
            prod = 1
            for r in 0:k
                prod *= 1/(l-r+im*ζ)
            end
            term = real(cis(ζ*χ*t)*(ψ+2(l+im*ζ)/ξ-ξ/(2(l+im*ζ+1)))*prod)
            #t1 = cos(ζ*(c1+c2*am)*t/m)*real((ψ+2(l+im*ζ)/ξ-ξ/(2(l+im*ζ+1)))*prod)
            #t2 = sin(ζ*(c1+c2*am)*t/m)*imag((ψ+2(l+im*ζ)/ξ-ξ/(2(l+im*ζ+1)))*prod)
            sum2 += ((-1)^l)*exp(2*l*χ*t)*term/(factorial(k-l)*factorial(l))
        end
        sum1 += exp(-2*k*χ*t)*((ξ/2)^(2*k+1))*sum2
    end
    return sum1
end

y1(t) = -(c1/(2*c2))*t+(m/c2)*log(y1Sum(t))

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

    tp = quadDragSim(dt = .0001)[5]
    global τ = tp + ep

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
