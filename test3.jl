using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
using QuadGK
include("projwDrag.jl")


function transform(f1)
    return function f2(x)
        return x[1] * f1(x[1] * x[2]) / f1(x[1])
    end
end

f1(t) = t

f2 = transform(f1)

println(f2([2,2]))

