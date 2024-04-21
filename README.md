# Marcus-Hush-Chidsey DOS rates in twisted Trilayer Graphene system

Implement MHC-DOS kinetics from the Julia-based [ElectrochemicalKinetics.jl](https://github.com/BattModels/ElectrochemicalKinetics.jl) package. Needs density of states of the solid as input. We employ the low-energy momentum space model to get DOS for the twisted trilayer graphene system.  

## Code descriptions

1. `angles.py` takes in DOS data files and outputs `q_dict.mat` with twist angle data. Will be input for `script.jl`

2. `script.jl` outputs `.mat` file with kox (kox_list), kred (kred_list), $\theta_{12}$ (q12_list) and $\theta_{23}$ (q23_list) variables for specified parameters $A$ , $\lambda$ and $\eta$

where $\lambda$ = reorganization energy (eV), $\eta$ = applied overpotential (V) and $A$ = proportionality constant for MHC-DOS theory. The output file has a format: k_data_ $A$ _ $\lambda$ _ $\eta$.mat

3. `eta_run_script.jl` runs `script.jl` at a range of eta values. 

The output mat file has a format: `k_data_ $A$ _ $\lambda$ _ $\eta$.mat`

4. `/sweep/` folder contains the .mat DOS files of the tTLG system at a range of theta12 and theta23.

5. `/Eo_var/` folder contains kinetic rates for a range of $\eta$ and $E_{o}$ (formal potential of redox couple wrt electrode).

For Ruthenium Hexamine 




`sweep_dos.m` uses output .mat file to generate colormaps of k[red/ox] for given [prefactor]_[lambda]_[eta]

`/Ruhex_data/` folder contains k_data_1.0_0.82_[eta].mat and colormaps(k[ox/red]_[eta].png) for analyzing MHCD rates of Ruhex with tTLG.






