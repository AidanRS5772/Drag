using PlotlyJS
using SpecialFunctions
using TickTock
include("projwDrag.jl")

function GammaAll(x)
    if isinteger(x)
        x = convert(Int128,x)
        if x <= 21
            return factorial(x-1)
        elseif x <= 0
            return Inf64
        else
            return factorial(big(x-1))
        end
    elseif typeof(x) == Float64
        if x > 0
            out  = gamma(x)
            if isinf(out)
                return exp(big(lng(x)))
            else
                return gamma(x)
            end
        else
            out = -π/(sin(π*x)*gamma(1-x))
            if out != 0.0
                return out
            else
                return -π*exp(-lng(1-x))/sin(π*x)
            end
        end
    elseif typeof(x) == ComplexF64
        out = gamma(x)
        if isinf(out) || out == 0.0
            return exp(big(lng(x)))
        else
            return gamma(x)
        end
    end
end
