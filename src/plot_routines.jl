function plot_dmrg_rg_flow!(p, color_lookup::Dict{Float64,String}, ms_lookup::Dict{Int,Real}) 
    Us = get_values_of_U(;path="../data/pbc/DMRG/")

    # move U clostest to Uc to beginning of array to bring curve to the front
    if 3.275 in Us  
        Us = [x for x in Us if x != 3.275]
        Us = [3.275;Us]
    end

	ζs = zeros(length(Us))
	ζserr = zeros(length(Us))

	for (b,U) in enumerate(reverse(Us))
        
	    Lmin = 16
        # get color for this U 
        col = color_lookup[U]
        # load array of L values available for U
        L_array_RG = get_values_of_L("../data/pbc/DMRG/",U;Lmin=Lmin)
        # Fit F to the linear regime to get the Luttinger parameter
        Ks_lin_RG = zeros(Float64,length(L_array_RG)) 
        Kerr_lin_RG = zeros(Float64,length(L_array_RG))
        for (i,L) in enumerate(L_array_RG) 
            ls, Fs = load_dmrg_pbc(L,U) 
            xs = lin_factor_pbc(ls,L)
            idx = window(L) 
            (Ks_lin_RG[i], _), Kerr_lin_RG[i] ,  _ = fit_fluc_lin_xs_pbc_err(xs, Fs, L, idx)  
        end
        # get pairs of L values to calculate ζ from each pair and to plot K(V) for each pair as points to the flow (use consecutive L values for pairs)
        stride = 1
        L_pairs = [(L_array_RG[i],L_array_RG[i+stride]) for i in 1:length(L_array_RG)-stride]
        K_pairs = [(Ks_lin_RG[i],Ks_lin_RG[i+stride]) for i in 1:length(Ks_lin_RG)-stride]
        K_err_pairs = [(Kerr_lin_RG[i],Kerr_lin_RG[i+stride]) for i in 1:length(Kerr_lin_RG)-stride]
        # loop over pairs (largest first)
        for (i, ((L1, L2), (K1,K2), (K1_err,K2_err))) in enumerate(zip(reverse(L_pairs),reverse(K_pairs),reverse(K_err_pairs))) 
            # get flow line for this pair
            ζ = get_zeta_from_Ks(L1,L2,K1,K2) 
            # estimate error based on error in K
            ζerr = max(abs(ζ-get_zeta_from_Ks(L1,L2,K1+K1_err,K2+K2_err) ) , abs(ζ-get_zeta_from_Ks(L1,L2,K1-K1_err,K2-K2_err) )) 
 
            # get flow curve for this ζ
            Ks, Vs = get_flow_curve(ζ; dk=0.00001)
            # from the flow curve locate the points K1 = K(L1) and K2 = K(L2)
            V1 = Vs[argmin(abs.(Ks.-K1))] 
            V2 = Vs[argmin(abs.(Ks.-K2))] 
            # get ζ and draw the solid flow curve for largest L pair
            if i == 1  
                ζs[b] = ζ
                ζserr[b] = ζerr 
                label  = latexstring(f"{U:5.3f}") 
                # handle different x-regions for ζ > 1 and ζ < 1
                if ζ > 1
                    iV0 = findfirst(Vs .== 0.0)
                    plot!(p, Ks[(Ks.>=2)][1:iV0], Vs[(Ks.>=2)][1:iV0]; color = col, label=nothing)  
                else
                    #if U > 3.285 || U < 3.265
                        plot!(p, Ks, Vs;color = col) # , label=label 
                    #end
                end   
                if U > 3.285 || U < 3.265
                    add_arrow!(p,Ks[(Ks.>=2)], Vs[(Ks.>=2)], 0.35*2/(4pi) ,col) 
                    add_arrow!(p,Ks[(Ks.>=2)], Vs[(Ks.>=2)], 0.28*2/(4pi) ,col)
                end
                # for the largest L pair also draw point K2, V2
                scatter!(p, [K2],[V2], xerr=[K2_err]  ;msc=col,lc=col,nice_points(col)...,lw=0.5,label=nothing, marker=:c,ms=2 )
                scatter!(p, [K2],[V2]  ; nice_points(col)..., markersize=ms_lookup[L2],lw=1 , label=nothing) # , label=nothing
            end 
            # draw point K1, V1   
            scatter!(p, [K1],[V1], xerr=[K1_err]  ;msc=col,lc=col,nice_points(col)...,lw=0.5,label=nothing, marker=:c,ms=2 )
            scatter!(p, [K1],[V1] ; nice_points(col)..., label=nothing, markersize=ms_lookup[L1],lw=1) 
        end
			  
	end

    Us = reverse(Us)
    idx = reverse(sortperm(Us))
    return Us[idx], ζs[idx], ζserr[idx]
end

function add_arrow!(p,x,y,y0,color; debug_print=false) 
    is_nan = isnan.(y)
    x = x[.!is_nan]
    y = y[.!is_nan]

    i_0 = argmin(abs.(y.-y0))

    debug_print && println(i_0)
    debug_print && println(y[i_0])
    debug_print && println(y0)
    debug_print && println(y)

    if i_0 > length(x)-1
        return
    end

    x_1 = x[i_0]
    x_2 = x[i_0+1]
    y_1 = y[i_0]
    y_2 = y[i_0+1]  
    # Calculate the differences
    dx = x_2 - x_1
    dy = y_2 - y_1 
    # Compute the angle in radians
    angle_rad = atan(dy*(4pi)/2, dx) 
    # Convert the angle to degrees 
    angle_deg = angle_rad * (180 / π)  
    debug_print && println(angle_deg)
    if angle_deg < 0
        angle_deg += 15 
    else
        angle_deg -= 10 
    end
    # Plot the arrow
    annotate!(p,x_2,y_2, Plots.text("➤", 5, :vcenter, :left, color, rotation = angle_deg ))
    #plot!(p, [x_1,x_2], [y_1,y_2], arrow=arrow=(0.01), color=color, linewidth=0.5, ms=0.5 )
end
 

function draw_legend!(p, Us, ζs, color_lookup::Dict{Float64,String}, ms_lookup::Dict{Int,Real}; lh=0.1, lU=0.42,lζ=0.75,pad=0.08, legendfontsize=10, top_margin=0.03, bottom_margin=0.0, highlight_U=[] )	
    # 1 |    . -- .2.0. 0.9| 
    #   |    .    .   .    |
    #   |    .    .   .    |
    # 0 |    .    .   .    |    
    #   0   lh    lU  lζ   1
    plot!(p; xlim=(0,1),ylim=(0,1),xticks=nothing,yticks=nothing, framestyle=nothing, axis=([], false))
    annotate!(p, (lU+lζ)/2, 1.0-top_margin, text(L"U/J", legendfontsize, :center, :center))
    annotate!(p, (lζ+1)/2, 1.0-top_margin, text(L"\zeta", legendfontsize, :center, :center))

    ys = reverse(range(bottom_margin, stop=1-top_margin, length=length(Us)+1))[2:end]

    for (U,ζ,y) in zip(Us,ζs,ys) 
        col = color_lookup[U]
        if U in highlight_U
            lw = 2 
            annotate!(p, lU+0.001, y+0.001, text(latexstring(f"{U:5.3f}"), legendfontsize, :left, :vcenter ))
            annotate!(p, lζ+0.001, y+0.001, text(latexstring(f"{ζ:5.4f}"), legendfontsize, :left, :vcenter )) 
        else 
            lw = 1 
        end
        annotate!(p, lU, y, text(latexstring(f"{U:5.3f}"), legendfontsize, :left, :vcenter ))
        annotate!(p, lζ, y, text(latexstring(f"{ζ:5.4f}"), legendfontsize, :left, :vcenter ))
        plot!(p, [lh,lU-pad], [y,y], color=col, lw=lw) 
    end
end 


function plot_vumps_rg_flow!(p, color_lookup::Dict{Float64,String}, ms_lookup::Dict{Int,Real} ) 
    U_vumps = get_values_of_U_vumps(;path="../data/pbc/VUMPS/")
    for U in U_vumps
        col = color_lookup[U]
        ls,xs,Fs = load_dmrg_vumps(U) 
        K, Kerr = fit_fluc_vumps_full_err(ls,Fs) 

        scatter!(p,[K],[0.00],  ;msc=col,lc=col,markercolor=col,lw=1,label=nothing, marker=:c,ms=4 )
		scatter!(p,[K],[0.00], xerr=[Kerr] ;msc=col,lc=col,markercolor=col,lw=0.5,label=nothing, marker=:c,ms=2 )
		scatter!(p,[K],[0.00] ;msc=col,lc=col,markercolor=:white,msw=0,label=nothing, marker=:c,ms=3,elw=0.0,lw=0) 
    end
end

function plot_qmc_rg_flow!(p, color_lookup::Dict{Float64,String}, ms_lookup::Dict{Int,Real} ) 
    U_qmc = get_values_of_U_qmc(;path="../data/pbc/QMC/") 
	for U in U_qmc 
		col = color_lookup[U]
        # load the qmc data
		L_qmc = get_values_of_L_qmc(U)  
		Ks_lin_RG = zeros(Float64,length(L_qmc)) 
		Kerr_lin_RG = zeros(Float64,length(L_qmc))
		
		for (i,L) in enumerate(L_qmc) 
			ls,Fs,δFs = load_dmrg_qmc(U, L) 
			xs = lin_factor_pbc(ls,L)
            # we cannot fully fit to the ℓ≈ L/2 region as the qmc statistical error is largest there
			idx = round(Int,L/8):round(Int,L/3)   
            # for the uncertainty include also the region up to ℓ=L/2
			mask_err = round(Int,L/3):round(Int,L/2)
			(Ks_lin_RG[i], _), Kerr_lin_RG[i] ,  _ = fit_fluc_lin_xs_pbc_err(xs, Fs, L, idx;mask_err=mask_err)  
		end
        # form the pairs from consecutive L values to solve for ζ for each pair
        stride = 1
		L_pairs = [(L_qmc[i],L_qmc[i+stride]) for i in 1:length(L_qmc)-stride]
		K_pairs = [(Ks_lin_RG[i],Ks_lin_RG[i+stride]) for i in 1:length(Ks_lin_RG)-stride]
		K_err_pairs = [(Kerr_lin_RG[i],Kerr_lin_RG[i+stride]) for i in 1:length(Kerr_lin_RG)-stride] 
				
		for (i, ((L1, L2), (K1,K2), (K1_err,K2_err))) in enumerate(zip(reverse(L_pairs),reverse(K_pairs),reverse(K_err_pairs))) 
			ζ = get_zeta_from_Ks(L1,L2,K1,K2)  
			Ks, Vs = get_flow_curve(ζ; dk=0.00001)
	 
			V1 = Vs[argmin(abs.(Ks.-K1))] 
			V2 = Vs[argmin(abs.(Ks.-K2))] 
			
			scatter!(p,[K1],[V1],  ;msc=col,lc=col,markercolor=col,lw=1,label=nothing, marker=:x,markersize=ms_lookup[L1]*3/4  )
			scatter!(p,[K1],[V1], xerr=[K1_err] ;msc=col,lc=col,markercolor=col,lw=0.5,label=nothing, marker=:x,ms=2 )
			if i == 1
				scatter!(p,[K2],[V2],  ;msc=col,lc=col,markercolor=col,lw=1,label=nothing, marker=:x,markersize=ms_lookup[L2]*3/4 )
				scatter!(p,[K2],[V2], xerr=[K2_err] ;msc=col,lc=col,markercolor=col,lw=0.5,label=nothing, marker=:x,ms=2 )
			end
		end
		
	end	
end

function plot_obc_rg_flow!(p, color_lookup::Dict{Float64,String}, ms_lookup::Dict{Int,Real} ) 
    factor = 4.1
    window_fun_obc = L-> round(Int,L/(2*factor)-L/16):round(Int,L/(2*factor)) 
    error_fun_obc = (xs)->round(Int,2length(xs)÷(8*factor)):round(Int,2length(xs)÷(2factor))

    U_obc = get_values_of_U_obc(;path="../data/obc/DMRG") 
    println(U_obc)  
    for U in U_obc
        col = color_lookup[U] 

        L_array = get_values_of_L_obc("../data/obc/DMRG",U;Lmin=300)[1:2:end]  
        Ks_lin_RG = zeros(Float64,length(L_array)) 
		Kerr_lin_RG = zeros(Float64,length(L_array))
 
        for (i,L) in enumerate(L_array) 
			ls, Fs = load_dmrg_obc(L,U);
			xs = lin_factor_pbc(ls, L);  
            
            idx = window_fun_obc(L)
               
			(Ks_lin_RG[i], _), Kerr_lin_RG[i] ,  _ = fit_fluc_lin_xs_pbc_err(xs, Fs, L, idx, mask_err=error_fun_obc(xs))   
		end  
        # form the pairs from consecutive L values to solve for ζ for each pair
        stride = 2
        L_pairs = [(L_array[i],L_array[i+stride]) for i in 1:length(L_array)-stride]
        K_pairs = [(Ks_lin_RG[i],Ks_lin_RG[i+stride]) for i in 1:length(Ks_lin_RG)-stride]
        K_err_pairs = [(Kerr_lin_RG[i],Kerr_lin_RG[i+stride]) for i in 1:length(Kerr_lin_RG)-stride] 
                
        for (i, ((L1, L2), (K1,K2), (K1_err,K2_err))) in enumerate(zip(reverse(L_pairs),reverse(K_pairs),reverse(K_err_pairs))) 
            ζ = get_zeta_from_Ks(L1,L2,K1,K2)  
            Ks, Vs = get_flow_curve(ζ; dk=0.00001)

            V1 = Vs[argmin(abs.(Ks.-K1))] 
            V2 = Vs[argmin(abs.(Ks.-K2))] 
            
            
            scatter!(p,[K1],[V1], xerr=[K1_err] ;msc=col,lc=col,markercolor=col,lw=0.5,label=nothing, marker=:square,ms=1 )
            scatter!(p,[K1],[V1],  ;msc=col,lc=col,markercolor=col,lw=1,label=nothing, marker=:square,markersize=2 ,nice_points(col)... )
            if i == 1 
                scatter!(p,[K2],[V2], xerr=[K2_err] ;msc=col,lc=col,markercolor=col,lw=0.5,label=nothing, marker=:square,ms=1 )
                scatter!(p,[K2],[V2],  ;msc=col,lc=col,markercolor=col,lw=1,label=nothing, marker=:square,markersize=2 ,nice_points(col)...)
            end
        end
    end
end

function parametrice_cicle(n, r, θ_range = (0,2π))
    θ = range(θ_range..., length=n+1)[1:end-1]  
    x = r * cos.(θ)   
    y = r * sin.(θ)  
    return x, y
end
function plot_circle!(p, n; r=1.0, kwargs...)   
    highlight_color = "#FF6961"

    x,y = parametrice_cicle(1000, r)
    plot!(p,(x .-r/2)./2,y; ls=:dash, lw=1, color=:black)
    x,y = parametrice_cicle(n, r)
    scatter!(p,(x .-r/2)./2, y;marker=:circle, kwargs...)

    scatter!(p,(x[1:4] .-r/2)./2, y[1:4];marker=:circle, mswidth  =1, nice_points(highlight_color)...)

    annotate!(p,-0.2,1.35, Plots.text("➤", 7, :vcenter, :left, color, rotation = 180 ))
    annotate!(p,0.355,0.15, Plots.text("➤", 7, :vcenter, :left, color, rotation = -90 ))
end 

function draw_periodic_lattice!(p; n=13, r=1.0, m=3, kwargs...)
    
    #plot_circular_sector!(p, m, n, r=r)
    plot_circle!(p, n; r=r, kwargs...)

    x,y = parametrice_cicle(10000, r+0.36, (-0.01,π/2-0.1)) 
    plot!(p, (x .-r/2 .-0.15)./2   , y;   color=:black  )
    annotate!(0.3,1.1,text(L"\ell",11,:black))
end 

function draw_open_lattice!(p; n=13, l = 10.0, kwargs...)
    highlight_color = "#FF6961"

    x = 0:l/n:l
     
    plot!(x,repeat([0],length(x)); ls=:dash, lw=1, color=:black)
    scatter!(x,repeat([0],length(x));marker=:circle, kwargs...)
    scatter!(x[5:10],repeat([0],length(x[5:10]));marker=:circle , mswidth  =1, nice_points(highlight_color)...)

    y = 0.2
    plot!(p,[l/2-2.1,l/2+2.1],[y,y];  color=:black, lw=1  )
    annotate!(p,l/2,0.35,text(L"\ell",11,:black))
    annotate!(p,x[10],y+.003, Plots.text("➤", 7, :vcenter, :left, color, rotation = 0 ))
    annotate!(p,x[5],y-.006, Plots.text("➤", 7, :vcenter, :left, color, rotation = 180 ))
end

function draw_open_lattice_horizontal!(p; n=13, l = 10.0, kwargs...)
    highlight_color = "#FF6961"

    y = 0:l/n:l 
    x = -0.2
     
    plot!(repeat([0],length(y)),y; ls=:dash, lw=1, color=:black)
    scatter!(repeat([0],length(y)),y;marker=:circle, kwargs...)
    scatter!(repeat([0],length(y[1:6])), y[1:6];marker=:circle , mswidth  =1, nice_points(highlight_color)...)

    plot!(p,[x,x],[y[1]-0.1,y[6]+0.1];  color=:black, lw=1  )
    annotate!(p,-0.4,(y[1]+y[6])/2,text(L"\ell",11,:black))
    annotate!(p,x-.01,y[6], Plots.text("➤", 7, :vcenter, :left, color, rotation = 90 ))
    annotate!(p,x+.01,y[1], Plots.text("➤", 7, :vcenter, :left, color, rotation = -90 ))
end

function plot_K_supplement!(p,  Us, color_lookup; load_fun=load_dmrg_pbc,L_lookup=get_values_of_L, window_fun=window, mask_err = xs -> length(xs)÷4:length(xs))
    for U in Us  
        Ls = L_lookup(U)
        c = color_lookup[U]

        Ks_lin = zeros(Float64,length(Ls))
        Ks_lin_err = zeros(Float64,length(Ls))
        as_lin = zeros(Float64,length(Ls))
        for (i,L) in enumerate(Ls)
            ls, Fs = load_fun(L,U) 
            xs = lin_factor_pbc(ls,L)
            idx = window_fun(L)  
            (Ks_lin[i], as_lin[i]), Ks_lin_err[i], _ = fit_fluc_lin_xs_pbc_err(xs,Fs,L,idx,mask_err=mask_err(xs)) 
        end    
        
        scatter!(p, x_trafo(Ls), Ks_lin,yerr=Ks_lin_err; marker=:circle, msw=1, ms=5, nice_points(c)...) 
    end
end

function plot_fluctuations_supplement!(p, Us, L, color_lookup; load_fun=load_dmrg_pbc, first_fit_value=96÷16, legend=true, window_fun=window, first_point=1,mask_err=xs->length(xs)÷4:length(xs))
    for (i_c,U) in enumerate(Us)
        c = color_lookup[U] 
		ls, Fs = load_fun(L,U);
		xs = lin_factor_pbc(ls, L);  
		idx = window_fun(L)
		(Kfit, afit), Kerr ,  _ = fit_fluc_lin_xs_pbc_err(xs, Fs, L, idx)   
		 
		scatter!(p, xs[first_point:end],Fs[first_point:end]; marker=:circle, msw=1, ms=5, label=latexstring(f"\\ U={U:3.3f}"), nice_points(c)...)
		 
        plot!(p, xs[1:mask_err(xs)[1]], fluctuations_pbc_large_L_xs(xs[1:mask_err(xs)[1]], (Kfit, afit), (L,)),color=:black,linewidth=0.5,ls=:dash,label=nothing) 
        plot!(p, xs[mask_err(xs)[end]:end], fluctuations_pbc_large_L_xs(xs[mask_err(xs)[end]:end], (Kfit, afit), (L,)),color=:black,linewidth=0.5,ls=:dash,label=nothing) 
		plot!(p, xs[mask_err(xs)], fluctuations_pbc_large_L_xs(xs[mask_err(xs)], (Kfit, afit), (L,)),color=:black,linewidth=1,ls=:dash,label=nothing)    
		plot!(p, xs[idx], fluctuations_pbc_large_L_xs(xs[idx], (Kfit, afit), (L,)),color=:black,linewidth=2,label=nothing)   

		# legend:
        if legend
            scatter!(p,[-0.04],[0.9-0.05*(i_c-1)];ms=5,nice_points(c)...)
            annotate!(-0.04,0.9-0.05*(i_c-1),text(latexstring(f"$\\ \\ \\ {U:3.3f}$"),:left,10))
        end
    end
    legend && annotate!(-0.027,0.9+0.05 ,text(latexstring(f"$U/J$"),:left,10))
end

function plot_fitting_method_literature!(p, U, L; colors = get_colors())
    _, Fs = load_dmrg_obc_literature_comparison(L,U)
    ls = 1:L
    Fs = L%2==0 ? [Fs;reverse(Fs)[2:end]] : [Fs;reverse(Fs)]
    xs = lin_factor_obc(ls,L) 

    W =L÷8:L
    

    # fit left window r=4
    K_left4 = zeros(length(W))
    a_left4 = zeros(length(W))
    for (i,w) in enumerate(W) 
        if 4+w > L-1
            continue
        end
        window = 4:min(4+w,L-1) 
        sol,_ = fit_fluc_lin_xs_pbc(xs, Fs, L, window) 
        K_left4[i] = sol[1]
        a_left4[i] = sol[2]
    end
    # fit left window r=16 
    K_left16 = zeros(length(W))
    a_left16 = zeros(length(W))
    for (i,w) in enumerate(W) 
        if 16+w > L-1
            continue
        end
        window = 16:min(16+w,L-1)
        sol,_ = fit_fluc_lin_xs_pbc(xs, Fs, L, window) 
        K_left16[i] = sol[1]
        a_left16[i] = sol[2]
    end
    # fit center window 
    K_center = zeros(length(W))
    a_center = zeros(length(W))
    for (i,w) in enumerate(W) 
        if  (L÷2-ceil(Int,w/2)) < 1 || (L÷2+floor(Int,w/2)) > L-1
            continue
        end
        window = max(1,(L÷2-ceil(Int,w/2))):min(L-1,(L÷2+floor(Int,w/2))) 
        sol,_ = fit_fluc_lin_xs_pbc(xs, Fs, L, window) 
        K_center[i] = sol[1]
        a_center[i] = sol[2]
    end

    K_fit = (K_center[1] + K_left4[1])/2
    K_err = abs(K_center[1] - K_left4[1])/2

    a_fit = (a_center[1] + a_left4[1])/2

    step = 2
    
    plot!(p,[0,L],[K_fit,K_fit], ribbon=K_err , color=:lightgray) 

    scatter!(p,(L ./ W)[1:step:end], K_center[1:step:end];m=:circle,nice_points(colors[1])...)
    scatter!(p,(L ./ W)[1:step:end], K_left16[1:step:end];m=:diamond,nice_points(colors[2])...)
    scatter!(p,(L ./ W)[1:step:end], K_left4[1:step:end];m=:square,nice_points(colors[3])...)


    plot!(p,[0,L/W[1]],[K_fit,K_fit];ls=:dash,color=:black)#,legend=:bottomleft,legendfontsize=9,label=nothing,foreground_color_legend = nothing,background_color_legend=nothing)

    return K_fit, K_err, a_fit
end

function plot_fss_method_literature!(p, U, Ls_lit, Ls_more; colors = get_colors())
    Ls = [Ls_lit;Ls_more]
    iend_lit = length(Ls_lit)

    Ks = zeros(Float64,length(Ls))
    Ks_err = zeros(Float64,length(Ls))
    for (i,L) in enumerate(Ls)
        nmax = 6
        _, Fs = load_dmrg_obc_literature_comparison(L,U;nmax)
        ls = 1:L
        Fs = L%2==0 ? [Fs;reverse(Fs)[2:end]] : [Fs;reverse(Fs)]
        xs = lin_factor_obc(ls,L)

        w =L÷8 

        # fit left window r=4  
        window = 4:min(4+w,L-1) 
        sol,_ = fit_fluc_lin_xs_pbc(xs, Fs, L, window) 
        K_left4 = sol[1] 
        # fit left window r=16  
        window = 16:min(16+w,L-1)
        sol,_ = fit_fluc_lin_xs_pbc(xs, Fs, L, window) 
        K_left16 = sol[1] 
        # fit center window  
        window = max(1,(L÷2-ceil(Int,w/2))):min(L-1,(L÷2+floor(Int,w/2))) 
        sol,_ = fit_fluc_lin_xs_pbc(xs, Fs, L, window) 
        K_center = sol[1] 

        K_fit = (K_center + K_left4)/2
        K_err = abs(K_center - K_left4)/2

        Ks[i] = K_fit
        Ks_err[i] = K_err
    end

    a, Kinf, Kinf_err = py_fit_lin_err_fss(1.0./log.(Ls[1:iend_lit]), Ks[1:iend_lit], Ks_err[1:iend_lit])
    
    xs = 0.0:0.0001:maximum(1.0./log.(Ls[1:iend_lit])) 

    #scatter!(p,1.0./log.(Ls[iend_lit+1:end]), Ks[iend_lit+1:end], yerr=Ks_err[iend_lit+1:end];marker=:circle, msw=1, ms=5, nice_points(colors[1])...)
    scatter!(p,1.0./log.(Ls[iend_lit+1:end]), Ks[iend_lit+1:end];marker=:diamond, msw=1, ms=5, nice_points(colors[1])...)
 
    plot!(p, xs, a*xs .+ Kinf; ls=:dash, color=:black)

    scatter!(p,1.0./log.(Ls[1:iend_lit]), Ks[1:iend_lit], yerr=Ks_err[1:iend_lit];marker=:circle, msw=1, ms=5, nice_points(colors[2])...)

    scatter!(p,[0.0], [Kinf], yerr=[Kinf_err];marker=:circle, msw=1, ms=5, nice_points("#000000")...)
    
    println("Esimated K:" , Kinf, " ± ", Kinf_err)
end