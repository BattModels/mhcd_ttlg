using Pkg; Pkg.activate("/Users/mbabar/Desktop/PhD/Analysis/Gerischer/")

loc = "/Users/mbabar/Desktop/PhD/Analysis/Gerischer/ElectrochemicalKinetics.jl/src/"

# Load libs
using MAT
using CSV
using DelimitedFiles
using Glob
using Interpolations
using Statistics
using LinearAlgebra
using QuadGK
include(loc*"dos.jl")
export DOS
using .DOS: DOSData, get_dos
export DOSData, get_dos
include(loc*"quantum_capacitance.jl")
include(loc*"kinetic_models.jl")
export MarcusHushChidseyDOS

## Numerator of the derived uncertainty expression
function int_numer(
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
function int_denom(
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

function compute_integral(
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
ox = true # Oxidation/Reduction
Vq_max = 0.6
Vq_min = -0.6
C_dl = 10.0
lambda = 0.82
Eo = -0.07
A = 1.0
kT = 0.026 #eV

## Get uncertainty in rates
# For single twist set calculation
single_calc = true
η_list = [0.0]
η = η_list[1]
end_str = "_kcut_6.9282_qtype_1_nq_484_zip.mat";
xcoord, ycoord = 1.25, 1.25 # Choose twist angles
xid, yid = findmin(abs.(q12_list .- xcoord))[2], findmin(abs.(q23_list .- ycoord))[2]
theta12, theta23 = chop_str(string(q12_list[xid])), chop_str(string(q23_list[yid]))
file = matopen(dir*"dos_q12_"*string(theta12)*"_q23_"*string(theta23)*end_str)
dos_tot = read(file, "dos_tot");
Elist = read(file, "E_list");
dos = [transpose(Elist) transpose(dos_tot)]
mhcd = MarcusHushChidseyDOS(A, lambda, dos, Ef=0.0) # Assuming Ef=0 (centered)
# Compute error
n, d = compute_integral(η, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max)
kox_err = n/d
n, d = compute_integral(η, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max)
kred_err = n/d
print(kox_err,",",kred_err,"\n")


## Multiple twists and eta
#η_list = [-1.0, -0.8, -0.6, -0.4, -0.2, 0.2, 0.4, 0.6, 0.8, 1.0]
#single_calc = false

mhcd = MarcusHushChidseyDOS(A, lambda, dos)
n, d = compute_integral(η, mhcd, ox; Vq_min=Vq_min, Vq_max=Vq_max)

kox_err = zeros(Float64, size(q12_list));
kred_err = zeros(Float64, size(q12_list));

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
        # Compute error
        n, d = compute_integral(η, mhcd, true; Vq_min=Vq_min, Vq_max=Vq_max)
        kox_err[i] = n/d
        n, d = compute_integral(η, mhcd, false; Vq_min=Vq_min, Vq_max=Vq_max)
        kred_err[i] = n/d
        print(i,",",kox_err[i],",",kred_err[i],"\n")
    end

    #print(k)

    # Write rate data as .mat
    file = matopen("err_data_"*string(A)*"_"*string(lambda)*"_"*string(η)*".mat", "w")
    write(file, "kox_err", kox_err)
    write(file, "kred_err", kred_err)
    write(file, "q12_list", q12_list)
    write(file, "q23_list", q23_list)
    close(file)
end

## Analyze output uncertainty data files

outdir = "DOS_unc_files/"
flist = glob(outdir*"k_err_"*"*")
matfile = matopen(flist[6]) # Any file number
kox_err = read(matfile, "kox_err")
kred_err = read(matfile, "kred_err")
q12_list = read(matfile, "q12_list")
q23_list = read(matfile, "q23_list")

xcoord, ycoord = 1.1, 4.7 # Choose twist angles
arr = [q12_list q23_list] .- [xcoord ycoord]
id = findmin(norm.(eachrow(arr)))[2] # OR eachcol(arr)
print("kox_ox=", kox_err[id],", k_red=", kred_err[id], "\n")

