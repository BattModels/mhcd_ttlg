using Pkg; Pkg.activate("/ocean/projects/cts180021p/mbabar/PhD/Gerischer/")

loc = "/ocean/projects/cts180021p/mbabar/PhD/Gerischer/ElectrochemicalKinetics.jl/src/"

using Distributed
using SharedArrays
using ClusterManagers

addprocs(62; exeflags="--project=/ocean/projects/cts180021p/mbabar/PhD/Gerischer/Project.toml")
@everywhere push!(LOAD_PATH, $loc)

# Load libs
@everywhere begin
    using Pkg; Pkg.activate("/ocean/projects/cts180021p/mbabar/PhD/Gerischer/")
    using MAT
    using CSV
    using DelimitedFiles
    using Glob
    using Interpolations
    using QuadGK
    using ElectrochemicalKinetics
    #calculate_Vdl_interp = ElectrochemicalKinetics.calculate_Vdl_interp
end

#include(loc*"rate_constant.jl")
#export compute_k, compute_k_cq

## Load functions
# Numerator of the derived uncertainty expression
@everywhere function int_numer(
    mhcd::MarcusHushChidseyDOS,
    ox::Bool;
    kT::Real = 0.026,
    η::Real = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( E .+ mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( E .- mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    #print(η,"\n")
    E -> mhcd.A .* marcus_term(E .+ η) .* fd(E)
end

## Denominator of the derived uncertainty expression
@everywhere function int_denom(
    mhcd::MarcusHushChidseyDOS,
    ox::Bool;
    kT::Real = 0.026,
    η::Real = 0.0,
    V_q::Real = 0.0,
)
    function marcus_term(E)
        local exp_arg
        if ox
            exp_arg = -(( E .+ mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        else
            exp_arg = -(( E .- mhcd.λ) .^ 2) ./ (4 * mhcd.λ * kT)
        end
        exp.(exp_arg)
    end
    fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    #print(η,"\n")
    E -> mhcd.A .* ((mhcd.dos.interp_func.(E .+ V_q)).^ 1) .* marcus_term(E .+ η) .* fd(E)
end

@everywhere function compute_integral(
    η,
    model::MarcusHushChidseyDOS,
    ox::Bool;
    Eo = -0.07, # E_f,red (solvent) - E_f,vac (bilayer)
    C_dl = 10.0,
    Vq_min = -0.5,
    Vq_max = 0.5,
    kT = 0.026,
    E_min = model.dos.E_min,
    E_max = model.dos.E_max,
)
    V_dl_interp = calculate_Vdl_interp(model.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_t = (Eo + η)
    V_dl = V_dl_interp(V_t)

    V_q = V_t - V_dl
    #print(V_dl,", ",V_q, ", ",V_t, "\n")
    if V_q < 0
           E_max = E_max .- 0.05
           E_min = E_min .- V_q .+ 0.05
    elseif V_q > 0
           E_min = E_min .+ 0.05
           E_max = E_max  .- V_q .- 0.05
    end
    #print(η, ", ",V_q, ", ",E_min,", ", E_max,"\n")
    numer = quadgk(int_numer(model, ox; kT = kT, η = η), E_min, E_max)[1]
    denom = quadgk(int_denom(model, ox; kT = kT, η = η, V_q = V_q), E_min, E_max)[1]
    return numer, denom
end


# Modify string -- helper function
@everywhere function chop_str(str::String)
    while str[length(str)] == '0'
          str = chop(str)
    end
    if str[length(str)] == '.'
       str = chop(str)
    end
    return str
end

# Load angles
file = matopen("../q_dict.mat")
q12_list = read(file, "q12");
q23_list = read(file, "q23");
q12_list = convert(SharedArray,collect(q12_list))
q23_list = convert(SharedArray,collect(q23_list))

# Load matlab dos data
dir = "../sweep/";
#end_str = "_kcut_6.9282_qtype_1_nq_484_zip.mat";
#file = matopen(dir*"dos_q12_"*string(theta12)*"_q23_"*string(theta23)*end_str)
#dos_tot = read(file, "dos_tot");
#Elist = read(file, "E_list");

# Get rate k
#dos = [transpose(Elist) transpose(dos_tot)]
## Define model
Vq_max = 0.6
Vq_min = -0.6
C_dl = 10.0
Eo = -0.07
λ = 0.82
A = 1.0
kT = 0.026 #eV
#η_list = [-0.3, -0.25, -0.2, -0.15, -0.1, -0.05, 0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
η_list = LinRange(-0.25, 0.25, 11)
η_list = convert(SharedArray,collect(η_list))

kox_err = zeros(Float64, size(q12_list));
kred_err = zeros(Float64, size(q12_list));
kox_err = convert(SharedArray,kox_err)
kred_err = convert(SharedArray,kred_err)

outdir = ""
for η in η_list
    @time begin
    print("\n","η = ",η,"\n")
    @sync @distributed for i in 1:size(q12_list)[1]
    	print(i,",")
        local file, theta12, theta23 
        theta12 = chop_str(string(q12_list[i]))
        theta23 = chop_str(string(q23_list[i]))
        mat_file = glob(dir*"dos_q12_"*theta12*"_q23_"*theta23*"*")[1]
        print("DOS file: "*mat_file,"\n")
        file = matopen(mat_file)
        dos_tot = read(file, "dos_tot");
        Elist = read(file, "E_list");
        dos = [transpose(Elist) transpose(dos_tot)]
        mhcd = MarcusHushChidseyDOS(A, λ, dos, Ef=0.0) # Assuming Ef=0 (centered)
        # Compute error
        n, d = compute_integral(η, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max)
        kox_err[i] = n/d
        n, d = compute_integral(η, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max)
        kred_err[i] = n/d
        print(kox_err[i],",",kred_err[i],"\n")
    end
    end
    # Write rate data as .mat
    file = matopen(outdir*"k_err_"*string(A)*"_λ_"*string(λ)*"_η_"*string(round(η,digits=2))*".mat", "w")
    write(file, "kox_err", collect(kox_err))
    write(file, "kred_err", collect(kred_err))
    write(file, "q12_list", collect(q12_list))
    write(file, "q23_list", collect(q23_list))
    #write(file, "prefactor", factor)
    close(file)
end

##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


