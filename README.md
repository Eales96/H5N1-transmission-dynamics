# H5N1-transmission-dynamics
Modelling of H5N1 transmission dynamics on dairy farms

This work was funded by the Livestock Biosecurity Fund - Cattle Compensation fund.
This code supports the work done under the agreement CCF 25.19: "Preparing for potential H5N1 incursions into Victorian dairy cattle".

# Reproducing the analyses
This code reproduces the analyses in the technical report produced under milestone 2 of the grant agreement. To reproduce the analyses all 'simulation' scripts should be run, and then all 'figure' scripts:

### Example simulation trajectories
To reproduce the figure showing an example infection trajectory for each transmission regime run the R script 'figures_example_simulation_timeseries.R'.
This will run three simulations (one for each transmission regime) and save the outputs in the 'simulations' subdirectory. It will then create the figure in the technical report. If the simulations have been run previosuly this section of code can be skipped when rerunning.

### Pre-emptive cohorting simulations
To reproduce all analyses of pre-emptive cohorting:
1) Run the R script 'simulations_preemptive.R'. This will save multiple 'parquet' files in the 'simulations' subdirectory.
2) Run the R script 'figures_preemptive.R'. This will then reproduce all relevant figures in the technical report.

### Bulk milk sample testing simulations
To reproduce all analyses for bulk milk sample testing:
1) Run the R script 'simulations_testing.R'. This will save multiple 'parquet' files in the 'simulations' subdirectory.
2) Run the R script 'figures_testing.R'. This will then reproduce all relevant figures in the technical report.

### Reactive cohorting simulations
To reproduce all analyses of reactive cohorting:
1) Run the R script 'simulations_reactive.R'. This will save multiple 'parquet' files in the 'simulations' subdirectory.
2) Run the R script 'figures_reactive.R'. This will then reproduce all relevant figures in the technical report.


# Underlying code
All functions used in performing the analyses are in the 'R' subdirectory.
The main functions for running simulations is the R scipt 'simulation.R'.
The R scipt 'simulation.R' simulates an outbreak on a dairy farm allowing the user to set as input:
1) df_fields: The structure and intial conditions of the dairy farm. Please see example uses of the function 'initialise_fields' in the main simulation scripts for how this function can be used. Note that the function initialises susceptible cows and so infections need to be seeded manually.
2) df_milk: The structure and intial conditions of the milking units. Please see example uses of the function 'initialise_milking_machines' in the main simulation scripts for how this function can be used.
3) params: The parameters to be used in the transmission model.  Please see examples of the parameters being defined in the main simulation scripts.
4) milking_time: The time of day at which the first milking occurs. Automatically sets to 6/24 (i.e. 6am) unless manually set.
5) milking_frequency: The frequency at which cattle are milked. Note that the milking times are then defined to be at evenly spaced intervals (which may not be appropriate for higher frequnecies of milking). Automatically sets to 2 which implies milking times of 6/24 and 18/24.
6) cleaning_time: the milking cohort before which enhanced cleaning occurs. Automatically sets to NA which implies no enhanced cleaning. Note that if simulating cleaning the cleaning_effect has to be set in params.

The code simulates the outbreak and will continue until there are no cattle in the exposed or infectious compartments, and it has been five days since the last milking unit was contaminated. The code returns a data.frame describing the number of cattle in each compartment (for each field) immediately following each milking period (i.e. time=10.5 is first milking on tenth day and time=11 is second milking on tenth day). There is also the convenient function 'run_multi_sims.R' that can be used to run multiple simulations using the 'simulation.R' function.

Note that there are analagous functions for performing reactive cohorting (i.e. 'simulation_reactive.R' and 'run_multi_sims_reactive.R'). This allows a 'reaction_threshold' and 'N_fields' to be set; when the total cumulative infections reach the reaction threshold (or higher) following a milking period the farm will be split into the number of specified field (N_fields) randomly. There is also analgous functions for performing simulations while storing the timing of all infection events ('simulation_detailed.R' and 'run_mulit_sims_detailed.R') these are necessary for simulating testing strategies (the exact timing of infections are required for simulating Ct value trajectories). 



