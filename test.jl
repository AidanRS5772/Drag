using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs()

sim = quadDragSim(dt = 10^-6 , track = false)

time = sim[3]

x(t) = (2*m/(3*c2))*log(imag(conj(r)*Whittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t))))+(lambda_p/(2*phi_p))*(1-exp(-phi_p*t))+(d_p0/6+g/(2*phi_p))*t 
y(t) = (2*m/(3*c2))*log(imag(conj(r)*Whittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t))))-(lambda_p/(2*phi_p))*(1-exp(-phi_p*t))-(d_p0/6-g/(2*phi_p))*t


T = LinRange(0,time,1000)
xa = []
ya = []
for t in T
    push!(xa,x(t))
    push!(ya,y(t))
end

pax = scatter(x = T , y=xa , mode = "line" , name = "x")
pay = scatter(x = T , y=ya , mode = "line" , name = "y")

plot([pax , pay])