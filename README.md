# MHCD-tTLG : Marcus-Hush-Chidsey DOS rates in Twisted Trilayer Graphene system

Implement MHC-DOS kinetics from the Julia-based [ElectrochemicalKinetics.jl](https://github.com/BattModels/ElectrochemicalKinetics.jl) package. Needs density of states of the solid as input. We employ the low-energy momentum space model to get DOS for the twisted trilayer graphene system. See code descriptions for instructions for data reproducibility. 

## Contact 

Mohammad Babar : mdbabar@umich.edu 

## Code descriptions

1. `angles.py` takes in DOS data files and outputs `q_dict.mat` with twist angle data. Will be input for `script.jl`

2. `script.jl` main Julia calculation script that outputs `.mat` file with oxidation rates (kox_list), reduction rates (kred_list), $\theta_{12}$ (q12_list) and $\theta_{23}$ (q23_list) variables for specified parameters $A$ , $\lambda$ and $\eta$.

where $\lambda$ = reorganization energy (eV), $\eta$ = applied overpotential (V) and $A$ = proportionality constant for MHC-DOS theory. Other input parameters:

i. `C_dl` : EDL capacitance (F)

ii. `V_dl` : EDL voltage (V)

iii. `C_q` : Quantum capacitance (F)

iv. `V_q` : Quantum capacitance voltage (V)

v. `Vq_min / Vq_max` : Min/Max range of Quantum capacitance voltage for interpolation (Eq. 3 in paper)

vi. `kT` : Thermal energy to temperature setting (0.26 eV at 300 K)

vii. `ef` : Fermi energy of the electrode


The output file has a format: `k_data_{A}_λ_{}_η_{}.mat`. Run Julia files with this command:

```
> julia script.jl
```

3. `eta_run_script.jl` runs `script.jl` at a range of eta values. 

The output mat file has a format: `k_data_{A}_λ_{}_η_{}.mat`.


4. `/sweep/` folder contains the .mat DOS files of the tTLG system at a range of $\theta_{12}$ and $\theta_{23}$. 

5. `/ttlg_dos/` repository contains data and instructions to generate twisted trilayer graphene DOS at a range of $\theta_{12}$ and $\theta_{23}$ (`1-5` degrees).

6. `/Eo_var/` folder contains kinetic rates for a range of $\eta$ and $E_{o}$ (formal potential of redox couple wrt electrode).

Formal potential of Ruthenium Hexamine, `E = -0.25 V` vs. Ag/AgCl electrode and reorg. energy `λ=0.82 eV` [Ref](https://www.nature.com/articles/s41557-021-00865-1).

Formal potential of twisted graphene, `E = -0.18 V` vs. Ag/AgCl electrode [Ref](https://www.nature.com/articles/s41557-021-00865-1). 

Hence `Eo = -0.25 - (-0.18) = -0.07 V` is used for Ruthenium Hexamine. The kinetic rate files for Figure 2b are stored in `/Eo_var/_0.07/` folder. 

Data for Figure 4 in paper is in `/Eo_var/0.3/` at equilibrium `k_data_1.0_λ_0.82_η_0.0.mat`.

6. `/trilayer_stacked/` repository contains Bernal stacked (ABA) trilayer graphene DOS data and MHC-DOS kinetic rates. The ABA rates are used as reference for color maps in paper (Fig. 2a, Fig.4, SI Fig. 3). The rate value is specifed in `sweep_dos.m`. See `/trilayer_stacked/README.md` for instructions.

7. `sweep_dos.m` is an analysis script that uses the rate file `.mat` to generate colormap of k $_{red/ox}$ or DOS with twist angles $\theta_{12}$ and $\theta_{23}$.

Specify surface vector `v` in lines 38-42 to either `kox_list` for oxidation rates, `kred_list` for reduction rates or `dos_max` for maximum DOS values (Figure 2a) as shown below.

```
x = q12_list;
y = q23_list;
v = kox_list;
v = kred_list;
v = dos_max;
```

`k_aba = 2.632e-5` on line 44 of `sweep_dos.m` is the equilibrium rate constant `ko` for ABA (Bernal) stacked trilayer graphene. It is used as reference for other twist angle rates. For calculation of `k_aba` see repository [trilayer_stacked](https://github.com/mbabar09/trilayer_stacked).

8. `/SI_fig3_data/` folder contains kinetic rate files at `η = 0.0` and `Eo = -0.07 V` for two reorganization energies `λ = 0.2 eV` and `λ = 1.2 eV`. These files can be loaded into `sweep_dos.m` to plot SI figure 3 color map. 


## Steps to reproduce data

1. Generate DOS data by following instructions in `/ttlg_dos/` repository. Data is stored in `/sweep/` folder.

2. Calculate reference ABA trilayer graphene DOS and corresponding kinetic rates in `/trilayer_stacked/` repository. Use `k_aba` value in `sweep_dos.m` 

3. Use `script.jl` or `eta_run_script.jl` files to generate kinetic rates in bulk for all twist angles (parallelized). 

4. Input rate file and run `sweep_dos.m` to visualize rates in a color map with $\theta_{12}$ and $\theta_{23}$ as `x` and `y` axes respectively.