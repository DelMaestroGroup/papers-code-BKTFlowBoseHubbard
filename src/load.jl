"""PyCall wrapper for numpy.loadtxt."""
function load_txt(path;delimiter="None")
    py""" 
    import numpy as np 
    
    def load_txt(path,delimiter):  
        if delimiter == "None":
            data = np.loadtxt(path, encoding="utf-8")
        else:
            data = np.loadtxt(path, encoding="utf-8", delimiter=delimiter)
    
        return data
    """
    
    return py"load_txt"(path,delimiter) 
end 


"""Loads fluctuation data from perioic boundary conditions DMRG simulations.

Parameters:
    L::Int: The length of the chain.
    U::Float64: The interaction strength.
    nmax::Int: The maximum number of bosons per site.
    path::String: The path to the folder containing the data (default: "../data/pbc/DMRG").
    
Returns: 
    ls::Vector{Float64}: The values of the subsystem size ℓ.
    Fs::Vector{Float64}: The values of the fluctuations F(ℓ).
"""
function load_dmrg_pbc(L::Int,U::Float64;nmax::Int=6,path::String="../data/pbc/DMRG")
    ls = collect(1.0:L//2)  
    Fs = try 
        filepath = joinpath(path, @sprintf "sigma2_L%02d_N%02d_nmax%02d_t+1.000_V+0.000_Vp+0.000_Usta%+4.3f_Uend%+4.3f_Unum0001.dat" L L nmax U U)
        Fs =  load_txt(filepath)[2:end]
    catch 
        filepath = joinpath(path, @sprintf "sigma2_L%02d_N%02d_nmax%02d_t+1.000_V+0.000_Vp+0.000_Usta%+4.5f_Uend%+4.5f_Unum0001.dat" L L nmax U U)
        Fs =  load_txt(filepath)[2:end]
    end
    return ls, Fs
end

function load_dmrg_vumps(U::Float64;nmax::Int=6,path::String="../data/pbc/VUMPS/") 
    filepath = joinpath(path, @sprintf "sigma2_BH_test_U%4.4f_3.txt" U)
    Fs = load_txt(filepath) 
    ls = collect(1.0:length(Fs))  
    xs = log.(ls)/(π^2)
    return ls, xs, Fs 
end

function load_dmrg_qmc(U::Float64,L::Int;nmax::Int=6,path::String="../data/pbc/QMC/") 
    filepath = joinpath(path, @sprintf "sigma2_aggregated_L%d_U%4.4f.npy" L U)
    data = npzread(filepath)
    ls, Fs, δFs = data[:,1], data[:,2], data[:,3]
    return ls, Fs, δFs 
end

function load_GP(;path="../data/pbc/GP/GP_Data_fulltrain.txt")
    data = load_txt(path;delimiter=",")
	U = data[:,1]
	z = 1e-5 * data[:,2]
	dz = 1e-5 * data[:,3]
	ζ = zeros(length(z))
	dζ = zeros(length(z))
	for ((i,_U), z_, dz_) in zip(enumerate(U),z,dz)
		if z_ < 0 
			dζ[i] = dz_/2
			ζ[i] = 1
			continue 
		end
		if _U < 3.273
			ζ[i] = 1 + sqrt(z_)
		else
			ζ[i] = 1 - sqrt(z_)
		end
		dζ[i] = dz_/(2*abs(1-ζ[i]))
	end
    return U, ζ, dζ
end

function load_GP_z(;path="../data/pbc/GP/GP_Data_fulltrain.txt")
    data = load_txt(path;delimiter=",")
	U = data[:,1]
	z = 1e-5 * data[:,2]
	dz = 1e-5 * data[:,3]
	 
    return U, z, dz
end

"""Checks directory of DMRG results and for a given U returns the values of L available."""
function get_values_of_L(path::String,U::Float64;Lmin::Int=0)  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))]
    Ls = Int[]
    for file in files 
        L = parse(Int, split(file,"_")[2][2:end])
        _U = parse(Float64, split(file,"_")[8][5:end])  
        if abs(U-_U) < 1e-4 && !(L  in Ls) && L >= Lmin
            push!(Ls,L)  
        end
    end
    sort!(Ls)
    return Ls
end

function get_values_of_L_qmc(U::Float64;Lmin::Int=0, path::String="../data/pbc/QMC/")  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))]
    Ls = Int[]
    for file in files 
        L = parse(Int, split(file,"_")[3][2:end])
        _U = parse(Float64, split(file,"_")[4][2:end-4])  
        if abs(U-_U) < 1e-4 && !(L  in Ls) && L >= Lmin
            push!(Ls,L)  
        end
    end
    sort!(Ls)
    return Ls
end

"""Checks directory of DMRG results and returns the values of U available."""
function get_values_of_U(;path::String="../data/pbc/DMRG/")  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))] 
    Us = Float64[]
    for file in files  
        U = parse(Float64, split(file,"_")[8][5:end])  
        if !(U in Us)
            push!(Us,U)     
        end
    end
    sort!(Us)
    return Us
end


function get_values_of_U_vumps(;path::String="../data/pbc/VUMPS/")  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))] 
    Us = Float64[]
    for file in files  
        U = parse(Float64, split(file,"_")[4][2:end])  
        if !(U in Us)
            push!(Us,U)     
        end
    end
    sort!(Us)
    return Us
end 

function get_values_of_U_qmc(;path::String="../data/pbc/QMC/")  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))] 
    Us = Float64[]
    for file in files   
        U = parse(Float64, split(file,"_")[4][2:end-4])  
        if !(U in Us)
            push!(Us,U)     
        end
    end
    sort!(Us)
    return Us
end 

function get_values_of_U_obc(;path::String="../data/obc/DMRG/")  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))] 
    Us = Float64[] 
    for file in files  
        flush(stdout)
        U = parse(Float64, split(file,"_")[9][5:end])  
        if !(U in Us)
            push!(Us,U)     
        end
    end
    sort!(Us)
    return Us
end 



# Supplemental Material 
function load_dmrg_obc(L::Int,U::Float64;nmax::Int=6,path::String="../data/obc/DMRG")
    ls = collect(1.0:L//2)  
    filepath = joinpath(path, @sprintf "sigma2_centerLR_L%02d_N%02d_nmax%02d_t+1.000_V+0.000_Vp+0.000_Usta%+4.3f_Uend%+4.3f_Unum0001_obc.dat" L L nmax U U)
    Fs =  load_txt(filepath)[2:end]
    return ls, Fs
end
function load_dmrg_obc_literature_comparison(L::Int,U::Float64;nmax::Int=4,path::String="../data/obc/DMRG_reproduce_literature")
    ls = collect(1.0:L//2)  
    filepath = joinpath(path, @sprintf "sigma2_L%02d_N%02d_nmax%02d_t+1.000_V+0.000_Vp+0.000_Usta%+4.3f_Uend%+4.3f_Unum0001_obc.dat" L L nmax U U) 
    Fs =  load_txt(filepath)[2:end]
    return ls, Fs
end
function get_values_of_L_obc(path::String,U::Float64;Lmin::Int=0)  
    files = [x for x in readdir(path) if isfile(joinpath(path,x))]
    Ls = Int[]
    for file in files 
        L = parse(Int, split(file,"_")[3][2:end])
        _U = parse(Float64, split(file,"_")[9][5:end])  
        if abs(U-_U) < 1e-4 && !(L  in Ls) && L >= Lmin
            push!(Ls,L)  
        end
    end
    sort!(Ls)
    return Ls
end