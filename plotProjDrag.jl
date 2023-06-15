using PlotlyJS
using Interpolations
using Printf

include("projwDrag.jl")

function pProjP()
    sim = qDSim(dt = 10.0^(-5) , track = false)
    simX = sim[1]
    simY = sim[2]
    simT = sim[3]

    ap = qDAproxP()
    apX = ap[1]
    apY = ap[2]
    apT = ap[3]
    xt = ap[4]
    yt = ap[5]
    Tt = ap[6]

    ptb = scatter(x = xt , y = yt , mode = "markers" , name = "Transition Points Both")
    pty = scatter(x = Tt , y = yt , mode = "markers" , name = "Transition Points Y")
    ptx = scatter(x = Tt , y = xt , mode = "markers" , name = "Transition Points X")
    ps = scatter(x=simX , y=simY , mode = "line" , name = "Sim")
    pa = scatter(x = apX , y=apY , mode = "line" , name = "Aprox")
    psx = scatter(x = simT , y = simX , mode = "line" , name = "Sim X")
    psy = scatter(x = simT , y = simY , mode = "line" , name = "Sim Y")
    pax = scatter(x = apT , y = apX , mode = "line" , name = "Aprox X")
    pay = scatter(x = apT , y = apY, mode = "line" , name = "Aprox Y")

    return [ps , pa , ptb] , [psx , psy , pax , pay , pty , ptx]
end

function pProjV()
    sim = qDSim(dt = 10.0^(-5) , track = false)
    simX = sim[1]
    simY = sim[2]
    simT = sim[3]
    dts = sim[4]

    simDX = []
    simDY = []
    simDT = simT[2:end-1]
    for n in 2:length(simT)-1
        push!(simDX , (simX[n+1]-simX[n-1])/(2*dts))
        push!(simDY , (simY[n+1]-simY[n-1])/(2*dts))
    end

    ap = qDAproxV()
    apDX = ap[1]
    apDY = ap[2]
    apDT = ap[3]
    xt = ap[4]
    yt = ap[5]
    Tt = ap[6]

    ptb = scatter(x = xt , y = yt , mode = "markers" , name = "Transition Points Both")
    pty = scatter(x = Tt , y = yt , mode = "markers" , name = "Transition Points Y")
    ptx = scatter(x = Tt , y = xt , mode = "markers" , name = "Transition Points X")
    ps = scatter(x=simDX , y=simDY , mode = "line" , name = "Simulation Velocity")
    pa = scatter(x = apDX , y=apDY , mode = "line" , name = "Aproximation Velocity")
    psdx = scatter(x = simDT , y = simDX , mode = "line" , name = "Sim DX")
    psdy = scatter(x = simDT , y = simDY , mode = "line" , name = "Sim DY")
    padx = scatter(x = apDT , y = apDX , mode = "line" , name = "Aprox DX")
    pady = scatter(x = apDT , y = apDY, mode = "line" , name = "Aprox DY")

    return [ps , pa , ptb] , [psdx , psdy , padx , pady , pty , ptx]
end

function errorRV1(vi , vf , N)

    function newton_method(f::Function, f_prime::Function, x0, tolerance=1e-8)
        x = x0
        prev_x = x + 2*tolerance  # Initialize prev_x to a value different from x
        while abs(x - prev_x) > tolerance
            prev_x = x
            x -= f(x) / f_prime(x)

            println("Newton's estimate:",x)
        end
        return x
    end


    errorData = []
    V = LinRange(vi , vf , N)
    for vn in V
        
        println("Velocity: ",vn)

        q = .677269

        vx1(t) = chi_p*exp(-c1*t/(2*c2))/(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))
        vy1(t) = -(c1/(2*c2))+(m*xi_p*omega_p/(2*c2))*exp(-omega_p*t)*((imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta+1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta)))-(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta-1))/imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))
        x1(t) = numInt(vx1 , 0, t)
        y1(t) = -(c1/(2*c2))*t+(m/c2)*log(abs(imag(conj(k)*Bessel(xi_p*exp(-omega_p*t),im*zeta))))
        rxy1(t) = vx1(t)/vy1(t)-q
    
        temp_u2(t) = exp((c1/m+c2*v2_p/(6*m))*t)*((abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))^(2/3))
        du2(t) = v2_m/temp_u2(t)-(g/temp_u2(t))*numInt(temp_u2 , 0 , t)
        u2(t) = numInt(du2 , 0 , t)
        #du2(t) = lambda_p*exp(-phi_p*t)-g/phi_p
        dv2(t) = (v2_p/3)+(2*m*phi_p/(3*c2))*((2*kappa_p-eta_p*exp(-phi_p*t))*real(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))-imag(conj(r)*(1+2*im*(kappa_p+mu_p))*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p*exp(-phi_p*t))))/imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))
        #u2(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
        v2(t) = (v2_p/3)*t+(4*m/(3*c2))*log(abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))
        x2(t) = (v2(t-t1)-u2(t-t1))/2 + d2_x
        y2(t) = (v2(t-t1)+u2(t-t1))/2 + d2_y
        vx2(t) = (dv2(t-t1)-du2(t-t1))/2
        vy2(t) = (dv2(t-t1)+du2(t-t1))/2
        ryx2(t) = vy2(t)/vx2(t)-q
    
        temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*t)*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*t)))
        temp_vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
        vy3(t) = v3_y/temp_y3(t-t2)-(g/temp_y3(t-t2))*numInt(temp_y3 , 0 , t-t2)
        y3(t) = numInt(temp_vy3 , 0 , t-t2)+d3_y
        vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t-t2)))*real(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))
        #vy3(t) = gammay*exp(-omega_0*(t-t2))-g/omega_0
        x3(t) = (v3_x/2)*(t-t2)+(m/c2)*log(abs(imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))))+d3_x
        #y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t-t2)))-(g/omega_0)*(t-t2)+d3_y
        ryx3(t) = -vy3(t)/vx3(t)-q
    
        temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*((abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3))
        dv4(t) = v4_p/temp_v4(t)-(g/temp_v4(t))*numInt(temp_v4 , 0 , t)
        v4(t) = numInt(dv4 , 0 , t)
        #dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
        du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
        #v4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
        u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))
        x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
        y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
        vx4(t) = (dv4(t-t3)-du4(t-t3))/2
        vy4(t) = (dv4(t-t3)+du4(t-t3))/2
        rxy4(t) = -vx4(t)/vy4(t)-q
    
        vx5(t) = chi_m*exp(-c1*(t-t4)/(2*c2))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
        vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
        x5(t) = numInt(vx5 , t4 , t) + d5_x
        y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(abs(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega)))+d5_y

        function ProjApX(t)
            if t < t1
                return x1(t)
            elseif t1 <= t < t2
                return x2(t)
            elseif t2 <= t < t3
                return x3(t)
            elseif t3 <= t < t4
                return x4(t)
            else t4 <= t 
                return x5(t)
            end     
        end

        function ProjApY(t)
            if t < t1
                return y1(t)
            elseif t1 <= t < t2
                return y2(t)
            elseif t2 <= t < t3
                return y3(t)
            elseif t3 <= t < t4
                return y4(t)
            else t4 <= t 
                return y5(t)
            end     
        end

        function ProjApDY(t)
            if t < t1
                return vy1(t)
            elseif t1 <= t < t2
                return vy2(t)
            elseif t2 <= t < t3
                return vy3(t)
            elseif t3 <= t < t4
                return vy4(t)
            else t4 <= t 
                return vy5(t)
            end     
        end

        instInputs(velocity = vn , theta = .9)

        chi_p = v0*cospi(θ/2)
        omega_p = (c1 + c2 * v0 * sinpi(θ/2)) / m
        xi_p = (c2 * chi_p) / (sqrt(2) * omega_p * m)
        zeta = sqrt(4*g*c2*m - c1^2) / (2*m*omega_p)
        k = -(π/(omega_p*sinh(π*zeta)))*(Bessel(xi_p , im*zeta)*((c2*v0*sinpi(θ/2)/m)+(c1/(2*m))-im*zeta*omega_p)+xi_p*omega_p*Bessel(xi_p,im*zeta-1))
        
        t1 = forwardEval(rxy1 , 0 , .1)

        d2_x = x1(t1)
        d2_y = y1(t1)
        v2_p = vy1(t1)+vx1(t1)
        v2_m = vy1(t1)-vx1(t1)
        phi_p = (2*c1+c2*v2_p)/(2*m)
        lambda_p = v2_m + g/phi_p
        eta_p = (3*c2*lambda_p)/(2*m*phi_p)
        kappa_p = (3*c2*g)/(4*m*phi_p^2)
        mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
        r = -(1/(2*eta_p*mu_p))*(rWhittaker(im*kappa_p,im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p))
        
        t2 = forwardEval(ryx2 , t1 , .01)
        
        d3_x = x2(t2)
        d3_y = y2(t2)
        v3_x = vx2(t2)
        v3_y = vy2(t2)
        omega_0 = (c1+c2*v3_x)/m
        gammay = v3_y+g/omega_0
        xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
        delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
        epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
        p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))
        
        t3 = forwardEval(ryx3 , t2 , .01)
        
        d4_x = x3(t3)
        d4_y = y3(t3)
        v4_p = vy3(t3)+vx3(t3)
        v4_m = vy3(t3)-vx3(t3)
        phi_m = (2*c1-c2*v4_m)/(2*m)
        lambda_m = v4_p + g/phi_m
        eta_m = (3*c2*lambda_m)/(2*m*phi_m)
        kappa_m = (3*c2*g)/(4*m*phi_m^2)
        mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
        s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

        t4 = forwardEval(rxy4 , t3 , .01)
        
        d5_x = x4(t4)
        d5_y = y4(t4)
        v5_y = vy4(t4)
        chi_m = vx4(t4)
        omega_m = (c1-c2*v5_y)/m
        xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
        Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
        l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
        l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

        sim = qDSim(track = false)
        simx = sim[1]
        simt = sim[3]
        simR = simx[end]
        simTime = simt[end]
        
        println("Sim Time: ",simTime)

        apT = newton_method(ProjApY , ProjApDY , simTime)
        apR = ProjApX(apT)

        error = abs(simR-apR)/simR
        println("Sim R: ",simR)
        println("Ap R: ",apR)
        push!(errorData , error)
    end

    return errorData
end

function errorRV2(vi , vf , N)

    function newton_method(f::Function, f_prime::Function, x0, tolerance=1e-8)
        x = x0
        prev_x = x + 2*tolerance  # Initialize prev_x to a value different from x
        while abs(x - prev_x) > tolerance
            prev_x = x
            x -= f(x) / f_prime(x)

            println("Newton's estimate:",x)
        end
        return x
    end


    errorData = []
    V = LinRange(vi , vf , N)
    for vn in V
        
        println("Velocity: ",vn)

        q = .677269

        temp_u2(t) = exp((c1/m+c2*v2_p/(6*m))*t)*((abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))^(2/3))
        du2(t) = v2_m/temp_u2(t)-(g/temp_u2(t))*numInt(temp_u2 , 0 , t)
        u2(t) = numInt(du2 , 0 , t)
        #du2(t) = lambda_p*exp(-phi_p*t)-g/phi_p
        dv2(t) = (v2_p/3)+(2*m*phi_p/(3*c2))*((2*kappa_p-eta_p*exp(-phi_p*t))*real(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))-imag(conj(r)*(1+2*im*(kappa_p+mu_p))*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p*exp(-phi_p*t))))/imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))
        #u2(t) = (lambda_p/phi_p)*(1-exp(-phi_p*t))-(g/phi_p)*t
        v2(t) = (v2_p/3)*t+(4*m/(3*c2))*log(abs(imag(conj(r)*exp(-im*mu_p*phi_p*t)*rWhittaker(im*kappa_p,im*mu_p,im*eta_p*exp(-phi_p*t)))))
        x2(t) = (v2(t)-u2(t))/2
        y2(t) = (v2(t)+u2(t))/2
        vx2(t) = (dv2(t)-du2(t))/2
        vy2(t) = (dv2(t)+du2(t))/2
        ryx2(t) = vy2(t)/vx2(t)-q
    
        temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*t)*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*t)))
        temp_vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
        vy3(t) = v3_y/temp_y3(t-t2)-(g/temp_y3(t-t2))*numInt(temp_y3 , 0 , t-t2)
        y3(t) = numInt(temp_vy3 , 0 , t-t2) + d3_y
        vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t-t2)))*real(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t-t2)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))
        #vy3(t) = gammay*exp(-omega_0*(t-t2))-g/omega_0
        x3(t) = (v3_x/2)*(t-t2)+(m/c2)*log(abs(imag(conj(p)*exp(-im*epsilon*omega_0*(t-t2))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t-t2))))))+d3_x
        #y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t-t2)))-(g/omega_0)*(t-t2)+d3_y
        ryx3(t) = -vy3(t)/vx3(t)-q
    
        temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*((abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3))
        dv4(t) = v4_p/temp_v4(t)-(g/temp_v4(t))*numInt(temp_v4 , 0 , t)
        v4(t) = numInt(dv4 , 0 , t)
        #dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
        du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
        #v4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
        u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))
        x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
        y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
        vx4(t) = (dv4(t-t3)-du4(t-t3))/2
        vy4(t) = (dv4(t-t3)+du4(t-t3))/2
        rxy4(t) = -vx4(t)/vy4(t)-q
    
        vx5(t) = chi_m*exp(-c1*(t-t4)/(2*c2))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
        vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
        x5(t) = numInt(vx5 , t4 , t) + d5_x
        y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(abs(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega)))+d5_y

        function ProjApX(t)
            if t < t2
                return x2(t)
            elseif t2 <= t < t3
                return x3(t)
            elseif t3 <= t < t4
                return x4(t)
            else t4 <= t 
                return x5(t)
            end     
        end

        function ProjApY(t)
            if t < t2
                return y2(t)
            elseif t2 <= t < t3
                return y3(t)
            elseif t3 <= t < t4
                return y4(t)
            else t4 <= t 
                return y5(t)
            end     
        end

        function ProjApDY(t)
            if t < t2
                return vy2(t)
            elseif t2 <= t < t3
                return vy3(t)
            elseif t3 <= t < t4
                return vy4(t)
            else t4 <= t 
                return vy5(t)
            end     
        end

        q = .677269

        instInputs(theta = .5 , velocity = vn)

        v2_p = v0*(sinpi(θ/2)+cospi(θ/2))
        v2_m = v0*(sinpi(θ/2)-cospi(θ/2))
        phi_p = (2*c1+c2*v2_p)/(2*m)
        lambda_p = v2_m + g/phi_p
        eta_p = (3*c2*lambda_p)/(2*m*phi_p)
        kappa_p = (3*c2*g)/(4*m*phi_p^2)
        mu_p = sqrt(12*c2*g*m*phi_p^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_p^2))/(4*m*phi_p^2)
        r = -(1/(2*eta_p*mu_p))*(rWhittaker(im*kappa_p,im*mu_p,im*eta_p)*((c2/(phi_p*m))*v2_p-im*(2*kappa_p-eta_p))+(1+im*2*(kappa_p+mu_p))*rWhittaker(im*kappa_p+1,im*mu_p,im*eta_p))
        
        t2 = forwardEval(ryx2 , 0 , .1)
        
        d3_x = x2(t2)
        d3_y = y2(t2)
        v3_x = vx2(t2)
        v3_y = vy2(t2)
        omega_0 = (c1+c2*v3_x)/m
        gammay = v3_y+g/omega_0
        xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
        delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
        epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
        p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))
        
        t3 = forwardEval(ryx3 , t2 , .1)
        
        d4_x = x3(t3)
        d4_y = y3(t3)
        v4_p = vy3(t3)+vx3(t3)
        v4_m = vy3(t3)-vx3(t3)
        phi_m = (2*c1-c2*v4_m)/(2*m)
        lambda_m = v4_p + g/phi_m
        eta_m = (3*c2*lambda_m)/(2*m*phi_m)
        kappa_m = (3*c2*g)/(4*m*phi_m^2)
        mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
        s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

        t4 = forwardEval(rxy4 , t3 , .1)
        
        d5_x = x4(t4)
        d5_y = y4(t4)
        v5_y = vy4(t4)
        chi_m = vx4(t4)
        omega_m = (c1-c2*v5_y)/m
        xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
        Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
        l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
        l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

        sim = qDSim(track = false)
        simx = sim[1]
        simt = sim[3]
        simR = simx[end]
        simTime = simt[end]
        
        println("Sim Time: ",simTime)

        apT = newton_method(ProjApY , ProjApDY , simTime)
        apR = ProjApX(apT)

        error = abs(simR-apR)/simR
        println("Sim R: ",simR)
        println("Ap R: ",apR)
        push!(errorData , error)
    end

    return errorData
end

function errorRV3(vi , vf , N)

    function newton_method(f::Function, f_prime::Function, x0, tolerance=1e-8)
        x = x0
        prev_x = x + 2*tolerance  # Initialize prev_x to a value different from x
        while abs(x - prev_x) > tolerance
            prev_x = x
            println("here 1")
            a = f(x)
            println("here 2")
            b = f_prime(x)
            println("here 3")
            x -= a/b

            println("Newton's estimate:",x)
        end
        return x
    end

    errorData = []
    V = LinRange(vi , vf , N)
    for vn in V
        
        println("Velocity: ",vn)

        q = .677269

        temp_y3(t) = exp((c1/m+c2*v3_x/(2*m))*t)*imag(conj(p)*exp(-im*epsilon*omega_0*t)*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*t)))
        temp_vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
        vy3(t) = v3_y/temp_y3(t)-(g/temp_y3(t))*numInt(temp_y3 , 0 , t)
        y3(t) = numInt(temp_vy3 , 0 , t)
        vx3(t) = (v3_x/2)+(m*omega_0/(2*c2))*((2*delta-xi_0*exp(-omega_0*(t)))*real(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t))))-imag(conj(p)*(1+2*im*(delta+epsilon))*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta+1,im*epsilon,im*xi_0*exp(-omega_0*(t)))))/imag(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t))))
        #vy3(t) = gammay*exp(-omega_0*(t))-g/omega_0
        x3(t) = (v3_x/2)*(t)+(m/c2)*log(abs(imag(conj(p)*exp(-im*epsilon*omega_0*(t))*rWhittaker(im*delta,im*epsilon,im*xi_0*exp(-omega_0*(t))))))
        #y3(t) = (gammay/omega_0)*(1-exp(-omega_0*(t)))-(g/omega_0)*(t)
        ryx3(t) = -vy3(t)/vx3(t)-q
    
        temp_v4(t) = exp((c1/m-c2*v4_m/(6*m))*t)*((abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))^(2/3))
        dv4(t) = v4_p/temp_v4(t)-(g/temp_v4(t))*numInt(temp_v4 , 0 , t)
        v4(t) = numInt(dv4 , 0 , t)
        #dv4(t) = lambda_m*exp(-phi_m*t)-g/phi_m
        du4(t) = (v4_m/3)-(2*m*phi_m/(3*c2))*((2*kappa_m-eta_m*exp(-phi_m*t))*real(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))-imag(conj(s)*(1+2*im*(kappa_m+mu_m))*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m*exp(-phi_m*t))))/imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))
        #4(t) = (lambda_m/phi_m)*(1-exp(-phi_m*t))-(g/phi_m)*t
        u4(t) = (v4_m/3)*t-(4*m/(3*c2))*log(abs(imag(conj(s)*exp(-im*mu_m*phi_m*t)*rWhittaker(im*kappa_m,im*mu_m,im*eta_m*exp(-phi_m*t)))))
        x4(t) = (v4(t-t3)-u4(t-t3))/2 + d4_x
        y4(t) = (v4(t-t3)+u4(t-t3))/2 + d4_y
        vx4(t) = (dv4(t-t3)-du4(t-t3))/2
        vy4(t) = (dv4(t-t3)+du4(t-t3))/2
        rxy4(t) = -vx4(t)/vy4(t)-q
    
        vx5(t) = chi_m*exp(-c1*(t-t4)/(2*c2))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
        vy5(t) = (c1/(2*c2))+(m*xi_m*omega_m/(2*c2))*exp(-omega_m*(t-t4))*(l_p*(Bessel(xi_m*exp(-omega_m*(t-t4)),Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),Omega+1))+l_m*(Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega-1)-Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega+1)))/(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega))
        x5(t) = numInt(vx5 , t4 , t) + d5_x
        y5(t) = (c1/(2*c2))*(t-t4)-(m/c2)*log(abs(l_p*Bessel(xi_m*exp(-omega_m*(t-t4)),Omega)+l_m*Bessel(xi_m*exp(-omega_m*(t-t4)),-Omega)))+d5_y

        function ProjApX(t)
            if t < t3
                return x3(t)
            elseif t3 <= t < t4
                return x4(t)
            else t4 <= t 
                return x5(t)
            end     
        end

        function ProjApY(t)
            if t < t3
                return y3(t)
            elseif t3 <= t < t4
                return y4(t)
            else t4 <= t 
                return y5(t)
            end     
        end

        function ProjApDY(t)
            if t < t3
                return vy3(t)
            elseif t3 <= t < t4
                return vy4(t)
            else t4 <= t 
                return vy5(t)
            end     
        end

        instInputs(theta = .3 , velocity = vn)
    
        v3_x = v0*cospi(θ/2)
        v3_y = v0*sinpi(θ/2)
        omega_0 = (c1+c2*v3_x)/m
        gammay = v3_y+g/omega_0
        xi_0 = (sqrt(2)*c2*gammay)/(m*omega_0)
        delta = (c2*g)/(sqrt(2)*m*(omega_0^2))
        epsilon = sqrt(2*(c2*g)^2-(c1*omega_0)^2)/(2*m*(omega_0^2))
        p = -(1/(2*xi_0*epsilon))*(rWhittaker(im*delta,im*epsilon,im*xi_0)*((c2*v3_x/(m*omega_0))-im*(2*delta-xi_0))+(1+2*im*(delta+epsilon))*rWhittaker(im*delta+1,im*epsilon,im*xi_0))
        
        t3 = forwardEval(ryx3 , 0 , .1)
        
        d4_x = x3(t3)
        d4_y = y3(t3)
        v4_p = vy3(t3)+vx3(t3)
        v4_m = vy3(t3)-vx3(t3)
        phi_m = (2*c1-c2*v4_m)/(2*m)
        lambda_m = v4_p + g/phi_m
        eta_m = (3*c2*lambda_m)/(2*m*phi_m)
        kappa_m = (3*c2*g)/(4*m*phi_m^2)
        mu_m = sqrt(-12*c2*g*m*phi_m^2+9*(c2^2)*(g^2)-4*(c1^2)*(phi_m^2))/(4*m*phi_m^2)
        s = (1/(2*eta_m*mu_m))*(rWhittaker(im*kappa_m,im*mu_m,im*eta_m)*((c2/(phi_m*m))*v4_m+im*(2*kappa_m-eta_m))-(1+im*2*(kappa_m+mu_m))*rWhittaker(im*kappa_m+1,im*mu_m,im*eta_m))

        t4 = forwardEval(rxy4 , t3 , .1)
        
        d5_x = x4(t4)
        d5_y = y4(t4)
        v5_y = vy4(t4)
        chi_m = vx4(t4)
        omega_m = (c1-c2*v5_y)/m
        xi_m = (c2*chi_m)/(sqrt(2)*omega_m*m)
        Omega = sqrt(4*g*c2*m+c1^2)/(2*m*omega_m)
        l_p = (π/(2*sinpi(Omega)))*(Bessel(xi_m,-Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))+Omega)+xi_m*Bessel(xi_m,-Omega-1))
        l_m = -(π/(2*sinpi(Omega)))*(Bessel(xi_m,Omega)*(((c1-2*c2*v5_y)/(2*m*omega_m))-Omega)+xi_m*Bessel(xi_m,Omega-1))

        sim = qDSim(track = false)
        simx = sim[1]
        simt = sim[3]
        simR = simx[end]
        simTime = simt[end]
        
        println("Sim Time: ",simTime)

        apT = newton_method(ProjApY , ProjApDY , simTime)
        apR = ProjApX(apT)

        error = abs(simR-apR)/simR
        println("Sim R: ",simR)
        println("Ap R: ",apR)
        push!(errorData , error)
    end

    return errorData
end