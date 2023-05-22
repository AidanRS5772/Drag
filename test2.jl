using PlotlyJS
using SpecialFunctions
using TickTock

function lng(z)
    m = 20
    b = Array{Float64}(undef, m+1)
    b[1] = 1.0
    b[2] = -0.5
    for n in 2:m
        if n%2 == 0
            local sum = 0
            for k in 0:n-1
                sum += (factorial(n)/(factorial(k)*factorial(n-k)))*(b[k+1]/(n-k+1))
            end
            b[n+1] = -sum
        else
            b[n+1] = 0
        end
    end
    b = b[3:2:m+1]

    sum = 0
    for k in 1:10
        sum += b[k]*z^(1-2*k)/((2*k)*(2*k-1))
    end
    return (z-1/2)*(log(abs(z))+im*angle(z))-z+(1/2)*log(2*π)+sum
end

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


function Gamma(x)
    if isinteger(x)
        x = convert(Int128,x)
        if x > 21
            return factorial(big(x-1))
        elseif x > 0
            return factorial(x-1)
        else
            return Inf
        end
    elseif typeof(x) == Float64
        if x > 0
            out = gamma(x)
            if isinf(out)||isnan(out)
                return exp(big(lng(x)))
            else
                return out
            end
        else
            out = π/(sinpi(x)*gamma(-x+1))
            if isinf(out)||isnan(out)
                return (-π/sinpi(x))*exp(-big(lng(-x+1)))
            else
                return out
            end
        end
    elseif typeof(x) == ComplexF64
        out = gamma(x)
        if isinf(out)|| abs(out) == 0
            return exp(big(lng(x)))
        else
            return out
        end
    end
end
