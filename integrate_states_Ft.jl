using Pkg; Pkg.activate("/Users/mbabar/Desktop/PhD/Analysis/Gerischer/")

loc = "/Users/mbabar/Desktop/PhD/Analysis/Gerischer/ElectrochemicalKinetics.jl/src/"

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

## Load functions
#Prefactor
function prefactor(
    mhcd::MarcusHushChidseyDOS;
    kT::Real = 0.026,
    V_q = 0.0,
)
    #fd(E) = ox ? 1 .- fermi_dirac(E; kT = kT) : fermi_dirac(E; kT = kT)
    Ft(E) = (1/4kT).*((sech.(E./2kT)).^2) # Thermal broadening function
    E -> ((mhcd.dos.interp_func.(E .+ V_q)).^ 1) .* Ft(E)
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

# Constants
Vq_max = 0.5
Vq_min = -0.5
C_dl=10.0
lambda = 0.82
A = 1.0
Vapp = 0.07
kT = 0.026 #eV
#mhcd = MarcusHushChidseyDOS(A, lambda, dos)

factor = zeros(Float64, size(q12_list));

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

    # Compute prefactor
    V_dl_interp = calculate_Vdl_interp(mhcd.dos.interp_func, Vq_min, Vq_max, C_dl)
    V_dl = V_dl_interp(Vapp)
    V_q = Vapp .- V_dl
    E_min = Elist[1,1];
    E_max = Elist[end,1];
    if V_q < 0
           E_min = E_min .- V_q
    elseif V_q > 0
           E_max = E_max .- V_q
    end
    #factor[i] = quadgk(prefactor(mhcd; kT = kT, V_q = V_q), E_min, E_max)[1]
    factor[i] = quadgk(prefactor(mhcd; kT = kT, V_q = V_q), -0.45, 0.45)[1]
    print(i,",",theta12,",",theta23,",",V_q,",",factor[i],"\n")

end

# Write prefactor to .mat
file = matopen("pref_"*string(kT)*"_"*string(Vapp)*".mat", "w")
write(file, "q12_list", q12_list)
write(file, "q23_list", q23_list)
write(file, "prefactor", factor)
close(file)


##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


