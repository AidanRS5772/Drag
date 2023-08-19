using PlotlyJS
using SpecialFunctions
using TickTock
using Cubature
using QuadGK
include("projwDrag.jl")



new_vx1(t) = chi_p*exp(-c1*t/(2*c2))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t) , im*zeta))
new_x1(t) = d1Int(new_vx1 , t)

temp_u2(t) = exp((c1/m+c2*v2_p/(6*m))*t)*(abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))^(2/3)
new_du2(t) = (v2_m-g*d1Int(temp_u2 , t))/temp_u2(t)
new_u2(t) = v2_m*d1Int(temp_u2 , t) - g*d2Int(triTorect(temp_u2) , (t , 1))

temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta , im*epsilon , im*xi_0*exp(-omega_0*t)))
new_vy3(t) = (v3_y-g*d1Int(temp_y3 , t))/temp_y3(t)
new_y3(t) = v3_y*d1Int(temp_y3 , t) - g*d2Int(triTorect(temp_y3) , (t , 1))

temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3)
new_dv4(t) = (v4_p - g*d1Int(temp_v4 , t))/temp_v4(t)
new_v4(t) = v4_p*d1Int(temp_v4 , t) - g*d2Int(triTorect(temp_v4) , (t , 1)) 

new_vx5(t) = chi_m*exp(-c1*t/(2*m))/(l_p*Bessel(xi_m*exp(-omega_m*(t)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t)),-Omega))
new_x5(t) = d1Int(new_vx5 , t)

instInputs()
preCalc(print = false)

oldx1 = []
oldvx1 = []
newx1 = []
newvx1 = []

oldu2 = []
olddu2 = []
newu2 = []
newdu2 = []

oldy3 = []
oldvy3 = []
newy3 = []
newvy3 = []

oldv4 = []
olddv4 = []
newv4= []
newdv4 = []

oldx5 = []
oldvx5 = []
newx5= []
newvx5 = []


T = LinRange(0 , 1 , 50)

#=
tick()
for t in T
    push!(oldx1 , x1(t))
    push!(oldvx1 , vx1(t))
end
tock()

tick()
for t in T
    push!(newx1 , new_x1(t))
    push!(newvx1 , new_vx1(t))
end
tock()
=#
tick()
for t in T
    push!(oldu2 , u2(t))
    push!(olddu2 , du2(t))
end
tock()

tick()
for t in T
    push!(newu2 , new_u2(t))
    push!(newdu2 , new_du2(t))
end
tock()
#=
tick()
for t in T
    push!(oldy3 , y3(t))
    push!(oldvy3 , vy3(t))
end
tock()

tick()
for t in T
    tick()
    push!(newy3 , new_y3(t))
    push!(newvy3 , new_vy3(t))
    tock()
end
tock()

tick()
for t in T
    push!(olddv4 , dv4(t))
    push!(oldv4 , v4(t))
end
tock()

tick()
for t in T
    tick()
    push!(newdv4 , new_dv4(t))
    push!(newv4 , new_v4(t))
    tock()
end
tock()

tick()
for t in T
    push!(oldx5 , x5(t))
    push!(oldvx5 , vx5(t))
end
tock()

tick()
for t in T
    push!(newx5 , new_x5(t))
    push!(newvx5 , new_vx5(t))
end
tock()
=#
#=
s1 = scatter(x = T , y = oldx1 , mode = "line" , name = "old x1")
s2 = scatter(x = T , y = oldvx1 , mode = "line" , name = "old vx1")
s3 = scatter(x = T , y = newx1 , mode = "line" , name = "new x1")
s4 = scatter(x = T , y = newvx1 , mode = "line" , name = "new vx1")
=#
s5 = scatter(x = T , y = oldu2 , mode = "line" , name = "old u2")
s6 = scatter(x = T , y = olddu2 , mode = "line" , name = "old du2")
s7 = scatter(x = T , y = newu2 , mode = "line" , name = "new u2")
s8 = scatter(x = T , y = newdu2 , mode = "line" , name = "new du2")
#=
s9 = scatter(x = T , y = oldy3 , mode = "line" , name = "old y3")
s10 = scatter(x = T , y = oldvy3 , mode = "line" , name = "old vy3")
s11 = scatter(x = T , y = newy3 , mode = "line" , name = "new y3")
s12 = scatter(x = T , y = newvy3 , mode = "line" , name = "new vy3")

s13 = scatter(x = T , y = oldv4 , mode = "line" , name = "old v4")
s14 = scatter(x = T , y = olddv4 , mode = "line" , name = "old dv4")
s15 = scatter(x = T , y = newv4 , mode = "line" , name = "new v4")
s16 = scatter(x = T , y = newdv4 , mode = "line" , name = "new dv4")

s17 = scatter(x = T , y = oldx5 , mode = "line", name = "old x5")
s18 = scatter(x = T , y = oldvx5 , mode = "line", name = "old vx5")
s19 = scatter(x = T , y = newx5 , mode = "line", name = "new x5")
s20 = scatter(x = T , y = newvx5 , mode = "line", name = "new vx5")
=#
#=
p1 = plot([s1 , s2 , s3 , s4])

p3 = plot([s9 , s10 , s11 , s12])
p4 = plot([s13 , s14 , s15 , s16])
p5 = plot([s17 , s18 , s19 , s20])
=#
p2 = plot([s5 , s6 , s7 , s8])
p2
