using PlotlyJS
using Printf
include("projwDrag.jl")

function newtonsWp(f::Function , g::Function , init , tol , preP::Bool)
    h = 10^(-12)
    px = init

    while tol < abs(g(px))
        px -= 2*g(px)*h/(g(px+h)-g(px-h))
    end

    if preP
        x = px - tol
    else
        x = px + tol
    end
    
    while tol < abs(f(x)-q)
        x -= 2*(f(x)-q)*h/(f(x+h)-f(x-h))
        if (x > px && preP) || (x < px && !preP)
            x = 2*px - x
        end
    end
    return x
end

function newtonNp(f::Function , init , tol)
    h = 10^(-12)
    x = init
    while tol < abs(f(x)-q)
        x -= 2*(f(x)-q)*h/(f(x+h)-f(x-h))
    end
    return x
end

instInputs(velocity = 20.0)
q = .677269

simdata = quadDragSim(dt = 10^(-5) , track = false)
xs = simdata[1]
ys = simdata[2]
time = simdata[3]
cnt = simdata[4]
dts = simdata[5]

TS = LinRange(0,time,cnt)
dxs = []
dys = []
for n in 2:cnt-1
    push!(dxs,(xs[n+1]-xs[n-1])/(2*dts))
    push!(dys,(ys[n+1]-ys[n-1])/(2*dts))
end

preCalc()
vx1(t) = chi_p*exp(-omega_p*t)
vy1(t) = -(c1/(2*c2))+(m*xi_p*omega_p/(2*c2))*exp(-omega_p*t)*((imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta+1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))-(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta-1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))



dta = 10^(-3)
TA = LinRange(0,time,floor(Int,time/dta))
xa = []
ya = []
dxa = []
dya = []
for t in TA
    push!(xa , x1(t))
    push!(ya , y1(t))
    push!(dxa , vx1(t))
    push!(dya , vy1(t))
end

pxs = scatter(x = TS , y = xs , mode = "line" , name = "X sim pos")
pys = scatter(x = TS , y = ys , mode = "line" , name = "Y sim pos")
pxa = scatter(x = TA , y = xa , mode = "line" , name = "X approx pos")
pya = scatter(x = TA , y = ya , mode = "line" , name = "Y approx pos")
pdxa = scatter(x = TA , y = dxa , mode = "line" , name = "X approx vel")
pdya = scatter(x = TA , y = dya , mode = "line" , name = "Y approx vel")
pdxs = scatter(x = TS[2:cnt-1] , y = dxs , mode = "line" , name = "X sim vel")
pdys = scatter(x = TS[2:cnt-1] , y = dys , mode = "line" , name = "Y sim vel")

rxs = abs.(dxs./dys)
rys = abs.(dys./dxs)
rxa = abs.(dxa./dya)
rya = abs.(dya./dxa) 

thresh = 5
rxs = map(x -> min(max(x,-thresh), thresh), rxs)
rxa = map(x -> min(max(x,-thresh), thresh), rxa)

pq = scatter(x = TS , y = fill(q, cnt) , mode = "line" , name = "q")
prxs = scatter(x = TS , y = rxs , mode = "line" , name = "X'/Y' sim")
prys = scatter(x = TS , y = rys , mode = "line" , name = "Y'/X' sim")
prxa = scatter(x = TA , y = rxa , mode = "line" , name = "X'/Y' approx")
prya = scatter(x = TA , y = rya , mode = "line" , name = "Y'/X' approx")


f1(t) = vx1(t)/vy1(t)
f2(t) = vy1(t)/vx1(t)
f3(t) = -vy1(t)/vx1(t)
f4(t) = -vx1(t)/vy1(t)

t1 = newtonsWp(f1 , vy1 , time/2 , dta^2 , true)
t2 = newtonNp(f2 , time/2 , dta^2)
t3 = newtonNp(f3 , time/2 , dta^2)
t4 = newtonsWp(f4 , vy1 , time/2 , dta^2 , false)

pt = scatter(x = [t1 , t2 , t3 , t4] , y = [f1(t1) , f2(t2) , f3(t3) , f4(t4)])

println(t1)
println(t2)
println(t3)
println(t4)

plot([pt , pq , prxa , prya])