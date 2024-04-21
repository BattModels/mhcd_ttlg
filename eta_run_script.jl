using Pkg; Pkg.activate("../ElectrochemicalKinetics.jl/")

loc = "../ElectrochemicalKinetics.jl/src/"

# Load libs
using MAT
using CSV
using DelimitedFiles
using Glob
using Interpolations
using QuadGK
include(loc*"dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos
include(loc*"quantum_capacitance.jl")
include(loc*"kinetic_models.jl")
export MarcusHushChidseyDOS
#include(loc*"rate_constant.jl")
#export compute_k, compute_k_cq

# Load functions
function integrand(
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

function compute_k_cq(
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
    k_rate = quadgk(integrand(model, ox; kT = kT, η = η, V_q = V_q), E_min, E_max)[1]
    return k_rate, V_q
end

# Modify string
function chop_str(str::String)
         while str[length(str)] == '0'
               str = chop(str)
         end
         if str[length(str)] == '.'
            str = chop(str)
         end
         return str
end

# Load angles
file = matopen("q_dict.mat")
q12_list = read(file, "q12");
q23_list = read(file, "q23");

# Load matlab dos data
dir = "sweep/";
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
lambda = 0.82
A = 1.0
η_list = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0]
kT = 0.026 #eV
#mhcd = MarcusHushChidseyDOS(A, lambda, dos)
#k_ox = compute_k_cq(0.4, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max)
#k_red = compute_k_cq(0.4, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max)

kox_data = zeros(Float64, size(q12_list));
kred_data = zeros(Float64, size(q12_list));
factor = zeros(Float64, size(q12_list));

for η in η_list
    print("\n","η = ",η,"\n")
    for i in 1:size(q12_list)[2]
        local file
        theta12 = chop_str(string(q12_list[i]))
        theta23 = chop_str(string(q23_list[i]))
        mat_file = glob(dir*"dos_q12_"*theta12*"_q23_"*theta23*"*")[1]
        file = matopen(mat_file)
        dos_tot = read(file, "dos_tot");
        Elist = read(file, "E_list");
        dos = [transpose(Elist) transpose(dos_tot)]
        mhcd = MarcusHushChidseyDOS(A, lambda, dos, Ef=0.0) # Assuming Ef=0 (centered)
        # Compute rates
        kox_data[i],V_q = compute_k_cq(η, mhcd, true; Eo=-0.07, Vq_min=Vq_min, Vq_max=Vq_max);
        kred_data[i],V_q = compute_k_cq(η, mhcd, false; Eo=-0.07, Vq_min=Vq_min, Vq_max=Vq_max);
        print(i,",",kox_data[i],",",kred_data[i],"\n")
    end

    #print(k)

    # Write rate data as .mat
    file = matopen("k_data_"*string(A)*"_"*string(lambda)*"_"*string(η)*".mat", "w")
    write(file, "kox_list", kox_data)
    write(file, "kred_list", kred_data)
    write(file, "q12_list", q12_list)
    write(file, "q23_list", q23_list)
    close(file)
end

##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


