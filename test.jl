using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")
include("plotProjDrag.jl")


instInputs()
#r = Whittaker(im*kappa_p,im*mu_p,im*eta_p)*((2*c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*Whittaker(im*kappa_p+1,im*mu_p,im*eta_p)


m = 20
b = Array{Float64}(undef, m+1)
b[1] = 1
b[2] = -1/2
for n in 1:m-1
    if n%2 == 1
        b[n+2] = 0
    else
        sum = 0
        for k in 0:n
            sum += factorial(n)*b[k+1]/(factorial(k)*factorial(n-k)*(n-k+2))
        end
        b[n+2] = sum*(n+1)
    end
end

println(b)
    
