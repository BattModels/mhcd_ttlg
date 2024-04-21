# Marcus-Hush-Chidsey DOS rates in twisted Trilayer Graphene system

Implement MHC-DOS kinetics from the Julia-based [ElectrochemicalKinetics.jl](https://github.com/BattModels/ElectrochemicalKinetics.jl) package. Needs density of states of the solid as input. We employ the low-energy momentum space model to get DOS for the twisted trilayer graphene system.  

## Code descriptions

`script.jl` # Outputs a .mat file with kox (kox_list), kred (kred_list), theta12 (q12_list) and theta23 (q23_list) variables
for given [prefactor]_[lambda]_[eta]

`eta_run_script.jl` runs script.jl at a range of eta values

# The output mat file has a format: k_data_<prefactor>_<lambda>_<eta>.mat
# where lambda = reorganization energy (eV), eta = applied overpotential (V)

`/sweep/` folder contains the .mat DOS files of	tTLG system at a range of theta12 and theta23.

`sweep_dos.m` uses output .mat file to generate colormaps of k<red/ox> for given <prefactor>_<lambda>_<eta>

`/Ruhex_data/` folder contains k_data_1.0_0.82_[eta].mat and colormaps(k<ox/red>_<eta>.png) for analyzing MHCD rates of Ruhex with tTLG.






