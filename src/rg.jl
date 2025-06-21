"""Integrand of RG flow equation."""
function invLogL_integrand!(F,K,ξ)  
    F .=  1.0 ./ K.^2 ./ (ξ .- log.(K/2.0) .- 2.0 ./ K)  
    nothing
end


function invLogL_integrand(ξ::Float64,K::Float64) 
    return  1.0 ./ K.^2 ./ (ξ-log.(K/2.0) .- 2.0 ./ K)
end 

"""Integral equation for the RG flow."""
function invLogLs(ξ, args)  
    L1,L2,K1,K2, prototype = args  
    prob = IntegralProblem{true}(BatchIntegralFunction(invLogL_integrand!, prototype), (K1,K2), (ξ))
    try 
        return solve(prob, QuadGKJL(order=50)).u + 2*(log(L1) - log(L2))
    catch 
        return Inf
    end
end

function invLogLs!(du, ξ, p)   
    du[1] = invLogLs(ξ[1], p) 
    nothing
end 


"""Solve the integral equation for ζ given K1, K2, L1, and L2."""
function get_zeta_from_Ks(L1,L2,K1,K2;ζstart = 0.7,ζstep=0.001, attempts=100, tol=1e-13, maxbisectionsteps = 1000000) 
    ζstart >= 1.0 && error("ζstart must be smaller than 1.0")

    args = (L1,L2,K1,K2, zeros(Float64,0) ) 
    ζfrom1 = 1.0 - ζstart
    # get first funtion value 
    fζ = invLogLs(ζstart,args)
    _fζ = invLogLs(ζstart + ζstep,args) 
    # check for initial monotonicity 
    sgn_m::Float64 = _fζ < fζ ? -1.0 : 1.0
    sgn::Float64 = 1.0
    # curve from -∞ is increasing and from +∞ is decreasing; decide on which site of the pole to start
    if sgn_m * fζ > 0.0 
        ζ = 1.0 + ζfrom1
        sgn = -1.0
        fζ = invLogLs(ζ,args)
        if fζ > 0.0
            # in this case, both sites are above, which may indicate that the solution is left to the starting point
            @warn "Could not determine starting point. Try to lower ζstart." 
            return NaN
        end
    else 
        # solution right of the starting point
        ζ = ζstart
        sgn = +1.0
    end  
    
    # counts failed attempts before failing 
    stopping_counter = 0 
    # find search region where zero is guaranteed to be between two values
    while true
        ζ += sgn * ζstep
        _fζ = invLogLs(ζ,args)
        # monotonicity check (jumping over the pole or numerical instability prevention)
        if sgn_m*sgn * _fζ < sgn_m*sgn * fζ
            # jumped over the pole or into numerical instability; attempt to jump back and increase stopping counter 
            ζ -= sgn * ζstep
            ζstep *= 0.5
            stopping_counter += 1
            if stopping_counter > attempts
                @warn "Could not find a zero. This may be due to numerical instability in the integrand for a zero close to the pole. You may try to incease the number of attempts or decrease the step size."
                return NaN
            end
            continue
        end
        # check if save region of zero is found
        if sgn_m*sgn * _fζ > 0.0 && sgn_m*sgn * fζ < 0.0
            break
        end 
        fζ = _fζ
    end
    # zero between ζ and ζ - sng * ζstep
    if _fζ < 0.0
        ζ1 = ζ
        ζ2 = ζ - sgn * ζstep
    else
        ζ1 = ζ - sgn * ζstep
        ζ2 = ζ
    end  
    # find zero using bisection method
    isteps = 0
    while abs(ζ2 - ζ1) > tol  
        ζ = 0.5 * (ζ1 + ζ2)
        fζ = invLogLs(ζ,args)
        if fζ > 0.0
            ζ2 = ζ
        else
            ζ1 = ζ
        end
        isteps += 1
        if isteps > maxbisectionsteps
            @warn "Bisection method did not converge to tol=$tol within $maxbisectionsteps steps. Try to increase the number of steps."
            return NaN
        end
    end 
    if abs(invLogLs(0.5 * (ζ1 + ζ2),args)) > 1e-8
        @debug "Residual: $(invLogLs(0.5 * (ζ1 + ζ2),args))"
        @warn "Bisection method did not converge to tol=1e-8. This may be due to a very steep function close to the zero."
        return NaN
    end
    return 0.5 * (ζ1 + ζ2) 
end

"""From a value of ζ, return a grid of K points with spacing dk and the corresponding V values according to the RG flow."""
function get_flow_curve(ξ::Float64 ; dk::Float64=0.001)
    K_array = collect(4.0:-dk:0.0) 
    # v/Er on y-axis 
    ys = sqrt(8)/π * sqrt.(0.0im .+ 2.0 ./ K_array .+  log.(K_array/2.0 .+ 0.0im) .-  ξ) 
    return K_array, real.(ys)   
end
 
 