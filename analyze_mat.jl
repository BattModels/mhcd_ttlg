using Pkg; Pkg.activate("/Users/mbabar/Desktop/PhD/Analysis/Gerischer/")

loc = "/Users/mbabar/Desktop/PhD/Analysis/Gerischer/ElectrochemicalKinetics.jl/src/"

# Load libs
using MAT
using CSV
using DelimitedFiles
using Glob
using LinearAlgebra
using Statistics

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

dos_dir = "sweep/";
file = matopen("q_dict.mat")
#end_str = "_kcut_6.9282_qtype_1_nq_484_zip.mat";
q12 = read(file, "q12");
q23 = read(file, "q23");

unc_dir = "DOS_unc_files/"
flist = glob(unc_dir*"k_err_"*"*")
matfile = matopen(flist[6]) # Any file number
kox_err = read(matfile, "kox_err")
kred_err = read(matfile, "kred_err")
q12_list = read(matfile, "q12_list")
q23_list = read(matfile, "q23_list")

#print("kox_ox=", kox_err[id],", k_red=", kred_err[id], "\n")

max_ox = 0
max_red = 0
for i in q12_list
    for j in q23_list
        xcoord, ycoord = i, j # Choose twist angles
        # if xcoord < 1.06 || ycoord < 1.06
        #     continue
        # end
        arr = [q12_list q23_list] .- [xcoord ycoord]
        id = findmin(norm.(eachrow(arr)))[2] # OR eachcol(arr)
        temp_ox, temp_red = kox_err[id], kred_err[id]

        arr2 = [transpose(q12) transpose(q23)] .- [xcoord ycoord]
        id2 = findmin(norm.(eachrow(arr2)))[2] # OR eachcol(arr)
        #xid, yid = findmin(abs.(q12 .- xcoord))[2], findmin(abs.(q23 .- ycoord))[2]
        theta12, theta23 = chop_str(string(q12[id2])), chop_str(string(q23[id2]))
        file = matopen(glob(dos_dir*"dos_q12_"*string(theta12)*"_q23_"*string(theta23)*"*")[1])
        #print(xcoord, ", ", ycoord, ". ", theta12, ", ", theta23, "\n")
        dos_tot = read(file, "dos_tot");
        Elist = read(file, "E_list");

        if temp_ox * mean(dos_tot)[1] > max_ox
            global max_ox = temp_ox * mean(dos_tot)[1]
            print(i, ", ", j, ", max_ox=", max_ox,"\n")
        end
    end
end