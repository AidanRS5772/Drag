using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")

instInputs()
out = eta_p*F11(im*mu_p-im*kappa_p+1/2,1+2*im*mu_p,im*eta_p)
println(log(imag(conj(r)*out)))
