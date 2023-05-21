using PlotlyJS
using SpecialFunctions
using TickTock




function Whittaker(a,b,x)
    sum = 0
    m = 20
    den = []
    num = []
    for n in 0:m
        push!(num , gamma(b-a+1/2+n))
        push!(den , exp(lng(2*b+1+n)))
        sum += gamma(b-a+1/2+n)*exp(-lng(2*b+1+n))*(x^n/factorial(n))
    end
    return exp(b*log(x))*sqrt(x)*exp(-x/2)*(exp(lng(2*b+1))/gamma(b-a+1/2))*sum , den , num
end

function oldWhittaker(a,b,x)
    sum = 0
    m = 20
    den = []
    num = []
    for n in 0:m
        push!(num , gamma(b-a+1/2+n))
        push!(den , gamma(2*b+1+n))
        sum += (gamma(b-a+1/2+n)/gamma(2*b+1+n))*(x^n/factorial(n))
    end
    return exp(-x/2)*(x^b)*sqrt(x)*(gamma(1+2*b)/gamma(b-a+1/2))*sum , den , num
end

n = 2090.64
m = 1045.55

