function fss_separatrix(Ls, (L,K))
    c = (1 + 2*log(L) - K *log(L))/(-2 + K)
    Ks = 1.0 ./ (c .+ log.(Ls)) .+ 2.0
    return Ks
end

function fss_superfluid(Ls, (L_points, K_points))
    # fit to the scaling law 
    function scaling_law(L, p,())
        #println(p)
        p[1] + p[2]*(1/L)^(2*(p[1]-2))
    end 
    (;sol, fit) = curve_fit(scaling_law, L_points, K_points, MVector{2}([2.5,0.5]), ())
    #println(sol)
    return [scaling_law(L, sol,()) for L in Ls]
end

function fss_mott_insulator(Ls, (L_points, K_points))
    # fit to the scaling law 
    function scaling_law(L, p,())
        #println(p)
        1 + p[1]*(1/L)^(2*(4-4*p[2]))
    end 
    (;sol, fit) = curve_fit(scaling_law, L_points, K_points, MVector{2}([0.5,0.99]), ())
    #println(sol)
    return [scaling_law(L, sol,()) for L in Ls]
end