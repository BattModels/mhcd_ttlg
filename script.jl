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
include(loc*"rate_constant.jl")
export compute_k, compute_k_cq

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
Vq_max = 0.5
Vq_min = -0.5
lambda = 0.82
A = 1.0
Vapp = -0.07
#mhcd = MarcusHushChidseyDOS(A, lambda, dos)
#k_ox = compute_k_cq(0.4, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max)
#k_red = compute_k_cq(0.4, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max)

kox_data = zeros(Float64, size(q12_list));
kred_data = zeros(Float64, size(q12_list));

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

for i in 1:size(q12_list)[2]
    local file
    theta12 = chop_str(string(q12_list[i]))
    theta23 = chop_str(string(q23_list[i]))
    mat_file = glob(dir*"dos_q12_"*theta12*"_q23_"*theta23*"*")[1]
    file = matopen(mat_file)
    dos_tot = read(file, "dos_tot");
    Elist = read(file, "E_list");
    dos = [transpose(Elist) transpose(dos_tot)]
    mhcd = MarcusHushChidseyDOS(A, lambda, dos)
    kox_data[i] = compute_k_cq(Vapp, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max);
    kred_data[i] = compute_k_cq(Vapp, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max);
    print(i,",",kox_data[i],",",kred_data[i],"\n")
end

#print(k)

# Write rate data as .mat
file = matopen("k_data_"*string(A)*"_"*string(lambda)*"_"*string(Vapp)*".mat", "w")
write(file, "kox_list", kox_data)
write(file, "kred_list", kred_data)
write(file, "q12_list", q12_list)
write(file, "q23_list", q23_list)
close(file)


##
#@time begin
#     mhcd = MarcusHushChidseyDOS(20.0, 0.82, dos)
#     k = compute_k(0.4, mhcd; calc_cq=true, Vq_min=-0.45, Vq_max=0.45)
#end


