using SpecialFunctions
using Cubature
using Printf

function dragEq(a, b)
    #put as the first entry the variable you want the equation to be repersenting, this 
    #does not include the inhomogeneous part of the y equation'
    return -(c1 / m) * a - (c2 / m) * a * sqrt(a^2 + b^2)
end

function qDSim(; dt=2.0^-16)
    time = 0
    t = []
    y = []
    x = []
    vY = []
    vX = []

    vy = v0 * sinpi(θ / 2)
    vx = v0 * cospi(θ / 2)

    ny = 0
    nx = 0

    while ny > -eps(Float64)

        time += dt
        push!(t, time)

        push!(x, nx)
        push!(y, ny)
        push!(vX, vx)
        push!(vY, vy)


        xk1 = dragEq(vx, vy)
        yk1 = dragEq(vy, vx) - g

        xk2 = dragEq(vx + dt * xk1 / 2, vy + dt * yk1 / 2)
        yk2 = dragEq(vy + dt * yk1 / 2, vx + dt * xk1 / 2) - g

        xk3 = dragEq(vx + dt * xk2 / 2, vy + dt * yk2 / 2)
        yk3 = dragEq(vy + dt * yk2 / 2, vx + dt * xk2 / 2) - g

        xk4 = dragEq(vx + dt * xk3, vy + dt * yk3)
        yk4 = dragEq(vy + dt * yk3, vx + dt * xk3) - g

        vx += (xk1 + 2 * xk2 + 2 * xk3 + xk4) * dt / 6
        vy += (yk1 + 2 * yk2 + 2 * yk3 + yk4) * dt / 6

        nx = x[end] + dt * vx
        ny = y[end] + dt * vy
    end

    ymaxi = findmax(y)[2]

    return x, y, t, vX, vY
end

function secantRF(f::Function, init ; track = false)
    if (track) println("\nSecant Root Finder:\n") end

    tol = 1e-8

    p0 = init
    p1 = init + 2 * tol
    p2 = 0

    fp0 = f(p0)
    fp1 = f(p1)
    if (track) println("t = " , init) end
    while abs(p1 - p0) > tol
        p2 = p1 - fp1 * ((p1 - p0) / (fp1 - fp0))
        if (track) println("t = " , p2) end
        p0 = p1
        p1 = p2

        fp0 = fp1
        fp1 = f(p2)
    end

    return p2
end

function d1Int(f::Function, t)
    out, _ = hquadrature(f, 0, t)
    return out
end

function d2Int(f::Function, (x1, x2))
    out, _ = hcubature(f, (0, 0), (x1, x2))
    return out
end

function triTorect(f1::Function)
    return function f2(x)
        return x[1] * f1(x[1] * x[2]) / f1(x[1])
    end
end

function recFunc(f::Function)
    return x -> 1 / f(x)
end

function besselCoef(z, x, n)
    if n <= 20
        out = (1 / (gamma(z + n + 1) * factorial(n))) * (x / 2)^(2 * n)
        if isinf(out) || isnan(out)
            return exp(-big(loggamma(z + n + 1))) * (big(x / 2)^n / factorial(n))
        else
            return out
        end
    else
        out = (1 / (gamma(z + n + 1) * factorial(big(n)))) * (x / 2)^(2 * n)
        if isinf(out) || isnan(out)
            return exp(-big(loggamma(z + n + 1))) * (big(x / 2)^n / factorial(big(n)))
        else
            return out
        end
    end
end

function Bessel(z, x)
    sum = 0
    n = 0

    while true
        val = ((-1)^n) * besselCoef(z, x, n)
        sum += val
        brp = abs(val)
        if brp < 1e-15
            break
        end
        n += 1
    end

    if isa(x, Complex) || isa(z, Complex)
        return convert(ComplexF64, (x / 2)^z * sum)
    else
        return convert(Float64, (x / 2)^z * sum)
    end
end

function rWhitCoef(a, b, x, n)
    if n <= 20
        out = (gamma(b - a + 1 / 2 + n) / gamma(2 * b + 1 + n)) * (x^n / factorial(n))
        if isinf(out) || isnan(out)
            return exp(big(loggamma(b - a + 1 / 2 + n) - loggamma(2 * b + 1 + n))) * (big(x)^n / factorial(n))
        else
            return out
        end
        return exp(loggamma(b - a + 1 / 2 + n) - loggamma(2 * b + 1 + n)) * (x^n / factorial(n))
    else
        out = (gamma(b - a + 1 / 2 + n) / gamma(2 * b + 1 + n)) * (x^n / factorial(big(n)))
        if isinf(out) || isnan(out)
            return exp(big(loggamma(b - a + 1 / 2 + n) - loggamma(2 * b + 1 + n))) * (big(x)^n / factorial(big(n)))
        else
            return out
        end
    end
end

function rWhittaker(a, b, x)
    sum = 0
    n = 0

    mult = exp(big(loggamma(2 * b + 1) - loggamma(b - a + 1 / 2)))

    while true
        val = rWhitCoef(a, b, x, n)
        sum += val
        brp = abs(val)
        if brp < (1e-15) / abs(mult)
            break
        end
        n += 1
    end

    return convert(ComplexF64, sqrt(x) * exp(-x / 2) * mult * sum)
end

function instInputs(; angle=0.9, velocity=10, diameter=0.2, mass=1)
    if angle >= 1 || angle <= 0
        DomainError("Launch angle must be in (0 , 1) and angle = $(angle)")
    end

    if v0 <= 0 
        DomainError("Lauch velocity must be greater then 0")
    end

    if D <= 0
        DomainError("Object size must be greater then 0")
    end

    if m <= 0
        DomainError("Object mass must be greater then 0")
    end

    global θ = convert(Float64 , angle)
    global v0 = convert(Float64 , velocity)
    global D = convert(Float64 , diameter)
    global m = convert(Float64 , mass)

    global c2 = (0.25) * D^2
    global c1 = (1.6e-4) * D
    global g = 9.8
    global q = 0.677269
end

function d1Vals(vx, vy; track=true)
    global χp = vx
    global ωp = (c1 + c2 * vy) / m
    global ξp = (c2 * χp) / (√2 * ωp * m)
    global ζ = √(4 * g * c2 * m - c1^2) / (2 * m * ωp)
    global k = -(π / (ωp * sinh(π * ζ))) * (Bessel(im * ζ, ξp) * (c2 * vy / m + c1 / (2 * m) - im * ζ * ωp) + ξp * ωp * Bessel(im * ζ - 1, ξp))

    if track
        println("χp = ", χp)
        println("ωp = ", ωp)
        println("ξp = ", ξp)
        println("ζ = ", ζ)
        println("k = ", k)
        println()
    end
end

vx1(t) = χp * exp(-(c1 / (2 * c2)) * t) / imag(conj(k) * Bessel(im * ζ, ξp * exp(-ωp * t)))
x1(t) = d1Int(vx1, t)
vy1(t) = -c1 / (2 * c2) - (m * ξp * ωp / c2) * exp(-ωp * t) * (imag(conj(k) * Bessel(im * ζ - 1, ξp * exp(-ωp * t))) - (ζ / ξp) * exp(ωp * t) * real(conj(k) * Bessel(im * ζ, ξp * exp(-ωp * t)))) / imag(conj(k) * Bessel(im * ζ, ξp * exp(-ωp * t)))
y1(t) = -(c1 / (2 * c2)) * t + (m / c2) * log(imag(conj(k) * Bessel(im * ζ, ξp * exp(-ωp * t))))
r1(t) = vx1(t) - q * vy1(t)

function d2Vals(x, y, vx, vy; track=true)
    global d2x = x
    global d2y = y
    global u2p = vy - vx
    global v2p = vy + vx
    global ϕp = (2 * c1 + c2 * v2p) / (2 * m)
    global λp = u2p + g / ϕp
    global ηp = 3 * c2 * λp / (2 * ϕp * m)
    global κp = 3 * c2 * g / (4 * m * ϕp^2)
    global μp = √(9 * c2^2 * g^2 + 12 * c2 * g * m * ϕp^2 - 4 * c1^2 * ϕp^2) / (4 * m * ϕp^2)
    global r = -(1 / (2 * ηp * μp)) * (rWhittaker(im * κp, im * μp, im * ηp) * (c2 * v2p / (ϕp * m) - im * (2 * κp - ηp)) + (1 + 2 * im * (μp + κp)) * rWhittaker(im * κp + 1, im * μp, im * ηp))

    if track
        println("d2x = ", d2x)
        println("d2y = ", d2y)
        println("u2p = ", u2p)
        println("v2p = ", v2p)
        println("ϕp = ", ϕp)
        println("λp = ", λp)
        println("ηp = ", ηp)
        println("κp = ", κp)
        println("μp = ", μp)
        println("r = ", r)
        println()
    end
end

temp_u2(t) = exp((c1 / m + c2 * v2p / (6 * m)) * t) * abs(imag(conj(r) * exp(-im * ϕp * μp * t) * rWhittaker(im * κp, im * μp, im * ηp * exp(-ϕp * t))))^(2 / 3)
du2(t) = (u2p - g * d1Int(temp_u2, t)) / temp_u2(t)
u2(t) = u2p * d1Int(recFunc(temp_u2), t) - g * d2Int(triTorect(temp_u2), (t, 1))

function dv2(t)
    M1 = conj(r) * exp(-im * ϕp * μp * t) * rWhittaker(im * κp, im * μp, im * ηp * exp(-ϕp * t))
    M2 = conj(r) * exp(-im * ϕp * μp * t) * rWhittaker(im * κp + 1, im * μp, im * ηp * exp(-ϕp * t))
    return v2p / 3 + (2 * m * ϕp / (3 * c2)) * ((2 * κp - ηp * exp(-ϕp * t)) * cot(angle(M1)) - imag((1 + 2 * im * (μp + κp)) * M2) / imag(M1))
end
v2(t) = v2p * t / 3 + (4 * m / (3 * c2)) * log(abs(imag(conj(r) * exp(-im * ϕp * μp * t) * rWhittaker(im * κp, im * μp, im * ηp * exp(-ϕp * t)))))


vx2(t) = (dv2(t - t1) - du2(t - t1)) / 2
vy2(t) = (dv2(t - t1) + du2(t - t1)) / 2
x2(t) = (v2(t - t1) - u2(t - t1)) / 2 + d2x
y2(t) = (v2(t - t1) + u2(t - t1)) / 2 + d2y
r2(t) = vy2(t) - q * vx2(t)

function d3Vals(x, y, vx, vy; track=true)
    global d3x = x
    global d3y = y
    global v3x = vx
    global v3y = vy
    global ω0 = (c1 + c2 * v3x) / m
    global γ = v3y + g / ω0
    global ξ0 = √2 * c2 * γ / (m * ω0)
    global δ = c2 * g / (√2 * m * ω0^2)
    global ϵ = √(2 * c2^2 * g^2 - c1^2 * ω0^2) / (2 * m * ω0^2)
    global p = -(1 / (2 * ξ0 * ϵ)) * (rWhittaker(im * δ, im * ϵ, im * ξ0) * (c2 * v3x / (m * ω0) - im * (2 * δ - ξ0)) + (1 + 2 * im * (δ + ϵ)) * rWhittaker(im * δ + 1, im * ϵ, im * ξ0))

    if track
        println("d3x = ", d3x)
        println("d3y = ", d3y)
        println("v3x = ", v3x)
        println("v3y = ", v3y)
        println("ω0 = ", ω0)
        println("γ = ", γ)
        println("ξ0 = ", ξ0)
        println("δ = ", δ)
        println("ϵ = ", ϵ)
        println("p = ", p)
        println()
    end
end

temp_y3(t) = exp((c1 / m + c2 * v3x / (2 * m)) * t) * abs(imag(conj(p) * exp(-im * ω0 * ϵ * t) * rWhittaker(im * δ, im * ϵ, im * ξ0 * exp(-ω0 * t))))
psvy3(t) = (v3y - g * d1Int(temp_y3, t)) / temp_y3(t)
psy3(t) = v3y * d1Int(recFunc(temp_y3), t) - g * d2Int(triTorect(temp_y3), (t, 1))

function psvx3(t)
    M1 = conj(p) * exp(-im * ω0 * ϵ * t) * rWhittaker(im * δ, im * ϵ, im * ξ0 * exp(-ω0 * t))
    M2 = conj(p) * exp(-im * ω0 * ϵ * t) * rWhittaker(im * δ + 1, im * ϵ, im * ξ0 * exp(-ω0 * t))
    return v3x / 2 + (m * ω0 / (2 * c2)) * ((2 * δ - ξ0 * exp(-ω0 * t)) * cot(angle(M1)) - imag((1 + 2 * im * (δ + ϵ)) * M2) / imag(M1))
end
psx3(t) = (v3x / 2) * t + (m / c2) * log(abs(imag(conj(p) * exp(-im * ω0 * ϵ * t) * rWhittaker(im * δ, im * ϵ, im * ξ0 * exp(-ω0 * t)))))

vx3(t) = psvx3(t - t2)
vy3(t) = psvy3(t - t2)
y3(t) = psy3(t - t2) + d3y
x3(t) = psx3(t - t2) + d3x
r3(t) = vy3(t) + q * vx3(t)

function d4Vals(x, y, vx; track=true)
    global d4x = x
    global d4y = y
    global u4p = -(1 + q) * vx
    global v4p = (1 - q) * vx
    global ϕm = (2 * c1 - c2 * u4p) / (2 * m)
    global λm = v4p + g / ϕm
    global ηm = 3 * c2 * λm / (2 * ϕm * m)
    global κm = 3 * c2 * g / (4 * m * ϕm^2)
    global μm = √(9 * c2^2 * g^2 - 12 * c2 * g * m * ϕm^2 - 4 * c1^2 * ϕm^2) / (4 * m * ϕm^2)
    global s = (1 / (2 * ηm * μm)) * (rWhittaker(im * κm, im * μm, im * ηm) * (c2 * u4p / (m * ϕm) + im * (2 * κm - ηm)) - (1 + 2 * im * (κm + μm)) * rWhittaker(im * κm + 1, im * μm, im * ηm))

    if track
        println("d4x = ", d4x)
        println("d4y = ", d4y)
        println("u4p = ", u4p)
        println("v4p = ", v4p)
        println("ϕm = ", ϕm)
        println("λm = ", λm)
        println("ηm = ", ηm)
        println("κm = ", κm)
        println("μm = ", μm)
        println("s = ", s)
        println()
    end
end

temp_v4(t) = exp((c1 / m - c2 * u4p / (6 * m)) * t) * imag(conj(s) * exp(-im * ϕm * μm * t) * rWhittaker(im * κm, im * μm, im * ηm * exp(-ϕm * t)))^(2 / 3)
dv4(t) = (v4p - g * d1Int(temp_v4, t)) / temp_v4(t)
v4(t) = v4p * d1Int(recFunc(temp_v4), t) - g * d2Int(triTorect(temp_v4), (t, 1))

function du4(t)
    M1 = conj(s) * exp(-im * ϕm * μm * t) * rWhittaker(im * κm, im * μm, im * ηm * exp(-ϕm * t))
    M2 = conj(s) * exp(-im * ϕm * μm * t) * rWhittaker(im * κm + 1, im * μm, im * ηm * exp(-ϕm * t))
    return u4p / 3 - (2 * m * ϕm / (3 * c2)) * ((2 * κm - ηm * exp(-ϕm * t)) * cot(angle(M1)) - imag((1 + 2 * im * (μm + κm)) * M2) / imag(M1))
end
u4(t) = u4p * t / 3 - (4 * m / (3 * c2)) * log(abs(imag(conj(s) * exp(-im * ϕm * μm * t) * rWhittaker(im * κm, im * μm, im * ηm * exp(-ϕm * t)))))

vx4(t) = (dv4(t - t3) - du4(t - t3)) / 2
vy4(t) = (dv4(t - t3) + du4(t - t3)) / 2
x4(t) = (v4(t - t3) - u4(t - t3)) / 2 + d4x
y4(t) = (v4(t - t3) + u4(t - t3)) / 2 + d4y
r4(t) = vx4(t) + q * vy4(t)

function d5Vals(x, y, vx ; track = true)
    global d5x = x
    global d5y = y
    global v5x = vx
    global v5y = -vx / q
    global ωm = (c1 - c2 * v5y) / m
    global ξm = c2 * v5x / (√2 * ωm * m)
    global Ω = √(4 * g * c2 * m + c1^2) / (2 * m * ωm)
    global lp = -(π / (2 * sin(π * Ω))) * (Bessel(-Ω, ξm) * ((c1 - 2 * c2 * v5y) / (2 * m * ωm) + Ω) + ξm * Bessel(-Ω - 1, ξm))
    global lm = (π / (2 * sin(π * Ω))) * (Bessel(Ω, ξm) * ((c1 - 2 * c2 * v5y) / (2 * m * ωm) - Ω) + ξm * Bessel(Ω - 1, ξm))

    if track 
        println("d5x = ", d5x)
        println("d5y = ", d5y)
        println("v5x = ", v5x)
        println("v5y = ", v5y)
        println("ωm = ", ωm)
        println("ξm = ", ξm)
        println("Ω = ", Ω)
        println("lp = ", lp)
        println("lm = ", lm)
        println()
    end
end

psvx5(t) = v5x * exp(-(c1 / (2 * c2)) * t) / (lp * Bessel(Ω, ξm * exp(-ωm * t)) + lm * Bessel(-Ω, ξm * exp(-ωm * t)))
psx5(t) = d1Int(psvx5, t)

psvy5(t) = (c1 / (2 * c2)) - (m * ξm * ωm / (2 * c2)) * exp(-ωm * t) * (lp * (Bessel(Ω + 1, ξm * exp(-ωm * t)) - Bessel(Ω - 1, ξm * exp(-ωm * t))) + lm * (Bessel(-Ω + 1, ξm * exp(-ωm * t)) - Bessel(-Ω - 1, ξm * exp(-ωm * t)))) / (lp * Bessel(Ω, ξm * exp(-ωm * t)) + lm * Bessel(-Ω, ξm * exp(-ωm * t)))
psy5(t) = (c1 * t / (2 * c2)) - (m / c2) * log(lp * Bessel(Ω, ξm * exp(-ωm * t)) + lm * Bessel(-Ω, ξm * exp(-ωm * t)))

vx5(t) = psvx5(t-t4)
vy5(t) = psvy5(t-t4)
x5(t) = psx5(t-t4) + d5x
y5(t) = psy5(t-t4) + d5y

function projP(t)
    if t <= t1
        return x1(t), y1(t)
    elseif t1 < t <= t2
        return x2(t), y2(t)
    elseif t2 < t <= t3
        return x3(t), y3(t)
    elseif t3 < t <= t4
        return x4(t), y4(t)
    elseif t4 < t
        return x5(t) , y5(t)
    end
end

function projV(t)
    if t <= t1
        return vx1(t), vy1(t)
    elseif t1 < t <= t2
        return vx2(t), vy2(t)
    elseif t2 < t <= t3
        return vx3(t), vy3(t)
    elseif t3 < t <= t4
        return vx4(t), vy4(t)
    else t4 < t
        return vx5(t), vy5(t)
    end
end

function projPx(t)
    if t <= t1
        return x1(t)
    elseif t1 < t <= t2
        return x2(t)
    elseif t2 < t <= t3
        return x3(t)
    elseif t3 < t <= t4
        return x4(t)
    elseif t4 < t
        return x5(t)
    end
end

function projPy(t)
    if t <= t1
        return y1(t)
    elseif t1 < t <= t2
        return y2(t)
    elseif t2 < t <= t3
        return y3(t)
    elseif t3 < t <= t4
        return y4(t)
    elseif t4 < t
        return y5(t)
    end
end

function projVx(t)
    if t <= t1
        return vx1(t)
    elseif t1 < t <= t2
        return vx2(t)
    elseif t2 < t <= t3
        return vx3(t)
    elseif t3 < t <= t4
        return vx4(t)
    else t4 < t
        return vx5(t)
    end
end

function projVy(t)
    if t <= t1
        return vy1(t)
    elseif t1 < t <= t2
        return vy2(t)
    elseif t2 < t <= t3
        return vy3(t)
    elseif t3 < t <= t4
        return vy4(t)
    else t4 < t
        return vy5(t)
    end
end

function preCalc(;print = true)
    if 1 > θ > 2*acot(q)/π
        d1Vals(v0 * cospi(θ / 2), v0 * sinpi(θ / 2) , track = print)
        global t1 = secantRF(r1, 0)
        if print println("t1 = " , t1 , "\n") end
        d2Vals(x1(t1), y1(t1), vx1(t1) , vy1(t1) , track = print)
        global t2 = secantRF(r2 , t1)
        if print println("t2 = " , t2 , "\n") end
        d3Vals(x2(t2) , y2(t2) , vx2(t2) , vy2(t2) , track = print)
        global t3 = secantRF(r3 , t2)
        if print println("t3 = " , t3 , "\n") end
        d4Vals(x3(t3) , y3(t3) , vx3(t3) , track = print)
        global t4 = secantRF(r4 , t3)
        if print println("t4 = " , t4 , "\n") end
        d5Vals(x4(t4) , y4(t4) , vx4(t4) , track = print)

    elseif 2*acot(q)/π >= θ > 2*atan(q)/π

        global t1 = 0

        d2Vals(0, 0 , v0 * cospi(θ / 2), v0 * sinpi(θ / 2) , track = print)
        global t2 = secantRF(r2 , 0)
        if print println("t2 = " , t2 , "\n") end
        d3Vals(x2(t2) , y2(t2) , vx2(t2) , vy2(t2) , track = print)
        global t3 = secantRF(r3 , t2)
        if print println("t3 = " , t3 , "\n") end
        d4Vals(x3(t3) , y3(t3) , vx3(t3) , track = print)
        global t4 = secantRF(r4 , t3)
        if print println("t4 = " , t4 , "\n") end
        d5Vals(x4(t4) , y4(t4) , vx4(t4) , track = print)

    elseif 2*atan(q)/π >= θ > 0

        global t1 = 0
        global t2 = 0

        d3Vals(0, 0 , v0 * cospi(θ / 2), v0 * sinpi(θ / 2) , track = print)
        global t3 = secantRF(r3 , 0)
        if print println("t3 = " , t3 , "\n") end
        d4Vals(x3(t3) , y3(t3) , vx3(t3) , track = print)
        global t4 = secantRF(r4 , t3)
        if print println("t4 = " , t4 , "\n") end
        d5Vals(x4(t4) , y4(t4) , vx4(t4) , track = print)
    end
end
    
