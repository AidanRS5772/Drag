using QuadGK
using SpecialFunctions
using Plots

x = LinRange(0,10,1000)
a = 1

y1 = []
iy2 = []
ry2 = []

for i in x
    push!(y1,1/gamma(i*im+a))
    integralr , _ = quadgk(t -> t^(-a)*exp(-t)*sin(log(t)*i), eps(Float64), 10, rtol=1e-16)
    rgr = ((-1)^a/π)*sinh(π*i)*integralr
    push!(ry2,rgr)
    integrali , _ = quadgk(t -> t^(-a)*exp(-t)*cos(log(t)*i), eps(Float64), 10, rtol=1e-16)
    rgi = ((-1)^a/π)*sinh(π*i)*integrali
    push!(iy2,rgi)
end

ry1 = real.(y1)
iy1 = imag.(y1)

p1 = plot(x,[ry1 iy1])
p2 = plot(x,[ry2 iy2])
display(plot(p1,p2))