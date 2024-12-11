"""Fit a model to data using NonlinearSolve.jl."""
function curve_fit(model, xdata, ydata, p0, params)
    data = (xdata, ydata, params)

    function lossfn!(du, u, data)
        (xs, ys, params) = data   
        du .= model.(xs, Ref(u), Ref(params)) .- ys
        return nothing
    end

    prob = NonlinearLeastSquaresProblem(
        NonlinearFunction(lossfn!, resid_prototype = similar(ydata)), p0, data)
    sol = solve(prob)
    u = sol.u
    fit = model.(xdata, Ref(u), Ref(params))
    return (;sol, fit)
end

"""PyCall based fitting routine using scipy.optimize.curve_fit. We use this to fit ζ(U) linearly near the critical point to obtain both U_c and the error in U_c."""
function py_fit_lin_err(U, ζs, ζserr)
    py""" 
    import numpy as np
    from scipy.optimize import curve_fit
    
    def fit_lin_err(x,y,yerr, p0):  
        popt,pcov = curve_fit(lambda U, Uc, a : a * (U - Uc) + 1.0 , x, y, method="lm",maxfev=100000,p0=p0,sigma=yerr,absolute_sigma=True)
        Uc, a = popt 
        Uc_err, a_err = np.sqrt(np.diag(pcov))

        popt,pcov = curve_fit(lambda U, invUc, a : a * (U - 1/invUc) + 1.0 , x, y, method="lm",maxfev=100000,p0=(1/Uc,a),sigma=yerr,absolute_sigma=True)
        invUc, inva = popt 
        invUc_err, inva_err = np.sqrt(np.diag(pcov))

        return Uc, Uc_err, invUc, invUc_err, a, a_err
        """

    Uc, dUc, invUc, dinvUc, α, dα = py"fit_lin_err"(U,ζs,ζserr, (3.3,1.0)) 
end

function py_fit_lin_err_fss(x, y, yerr)
    py""" 
    import numpy as np
    from scipy.optimize import curve_fit
    
    def fit_lin_err_fss(x,y,yerr, p0):  
        popt,pcov = curve_fit(lambda x, a, y0 : a*x+y0 , x, y, method="lm",maxfev=100000,p0=p0,sigma=yerr,absolute_sigma=True)
        a, y0 = popt 
        a_err, y0_err = np.sqrt(np.diag(pcov)) 
        return a, y0, y0_err
        """

        a, y0, y0_err = py"fit_lin_err_fss"(x,y,yerr, (0.1,2.0)) 
end

"""Window defining the region of largest ℓ≈L/2 that we fit the linear form to."""
function window(L)
    i_fit_start(L) = round(Int,L/2) - max(round(Int,L/16),2) 
	i_fit_end(L) = round(Int,L/2)
	return i_fit_start(L):i_fit_end(L)
end

"""Turn ℓ into x=1/(π^2)log(sin(πℓ/L)) for periodic boundary conditions to perfom linear fits to F(x)."""
function lin_factor_pbc(ls, L::Int)
    """
    F_pbc = lin_factor * K + A
    """
    return 1/(π^2) * log.(sin.(π*ls/L))
end

"""Turn ℓ into x=1/(2π^2)log(sin(πℓ/2L)) for periodic boundary conditions to perfom linear fits to F(x)."""
function lin_factor_obc(ls, L::Int)
    """
    F_obc = lin_factor * K + A
    """
    return 1/(2π^2) * real.(log.(Complex.(sin.(π*ls/L))))
end

"""Approximate LL form of the fluctuations for ℓ≈L/2 assuming large L."""
function fluctuations_pbc_large_L_xs(xs, params, p)
    K, A = params
    L, = p

    return K * xs .+ A 
end 

"""Fit fluctuations to obtain Luttinger parameter and error in K.

Parameters: 
    xs::Vector{Float64}: The transformed values of the subsystem size ℓ: \\frac{1}{2\\pi^2}\\ln \\sin\\frac{\\pi\\ell}{L}.
    Fs::Vector{Float64}: The values of the fluctuations F(ℓ).
    L::Int: The length of the chain.
    mask::UnitRange{Int}: The range of indices defining the ℓ values to fit the linear form to.
    mask_err::UnitRange{Int}: The range of ℓ values to estimate the error in K (error defined via the uncertainty in the fit procedure due to choice of fit interval).

Returns:
    sol::Vector{Float64}: The fitted parameters K and A.
    K_err::Float64: The estimate of the error in K.
    fit::Vector{Float64}: The fitted values of the fluctuations.
"""
function fit_fluc_lin_xs_pbc_err(xs, Fs, L, mask; mask_err = length(xs)÷4:length(xs))
    K, A = 2.0, 1.0 
    params = (L,)
    (;sol, fit) = curve_fit(fluctuations_pbc_large_L_xs, xs[mask], Fs[mask], MVector{2}([K,A]), params)
    
    K = sol[1]

    K_err = 0.0
    for i = 2:length(mask_err)-1
        _mask_err = mask_err[end-i:end]
        sol_err , _ = curve_fit(fluctuations_pbc_large_L_xs, xs[_mask_err], Fs[_mask_err], MVector{2}([K,A]), params)
        K_err = max(K_err, abs(K - sol_err[1]))
    end 
    return sol, K_err, fit
end 

function fit_fluc_lin_xs_pbc(xs, Fs, L, mask)
    K, A = 2.0, 1.0 
    params = (L,)
    (;sol, fit) = curve_fit(fluctuations_pbc_large_L_xs, xs[mask], Fs[mask], MVector{2}([K,A]), params) 
    return sol, fit
end 



"""Fluctuations LL L->∞ form for use with infinite VUMPS results."""
function fluc_td_limit_full(ℓ, K, a)
    return  K/(2*π^2) *  log.(1 + ℓ^2/a^2)  
end

function fluc_td_limit_full(ℓs::AbstractArray, K, a)
    return [fluc_td_limit_full(ℓ, K, a) for ℓ in ℓs]
end

function fluc_td_limit_full(ℓ, u0::AbstractArray, params::Tuple)
    K, a = u0
    return  fluc_td_limit_full(ℓ, K, a)
end 
 
function fit_fluc_vumps_full_err(ℓs,Fℓs; mask = (length(ℓs)-100):length(ℓs), mask_1=length(ℓs)-4:length(ℓs), mask_2 = 4:length(ℓs) )
    params = (  )
    # error via results for different fit window choices
    (K_1,a_1), _ = curve_fit(fluc_td_limit_full, ℓs[mask_1], Fℓs[mask_1], MVector{2}([2.0,1.0]), params) 
    (K_2,a_2), _ = curve_fit(fluc_td_limit_full, ℓs[mask_2], Fℓs[mask_2], MVector{2}([2.0,1.0]), params)
    
    (K,a), _ = curve_fit(fluc_td_limit_full, ℓs[mask], Fℓs[mask], MVector{2}([2.0,1.0]), params)
    return K, max(abs(K-K_1),abs(K-K_2)) 
end