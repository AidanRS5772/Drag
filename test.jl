using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

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

preCalc()
vx1(t) = chi_p*exp(-omega_p*t)
vy1(t) = -(c1/(2*c2))+(m*xi_p*omega_p/(2*c2))*exp(-omega_p*t)*((imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta+1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))-(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta-1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))



dxa = 10^(-3)
f1(t) = vx1(t)/vy1(t)
f2(t) = vy1(t)/vx1(t)
f3(t) = -vy1(t)/vx1(t)
f4(t) = -vx1(t)/vy1(t)

τ1 = newtonsWp(f1 , vy1 , time/2 , dxa^2 , true)
τ2 = newtonNp(f2 , time/2 , dxa^2)
τ3 = newtonNp(f3 , time/2 , dxa^2)
τ4 = newtonsWp(f4 , vy1 , time/2 , dxa^2 , false)

println(τ1)
println(τ2)
println(τ3)
println(τ4)