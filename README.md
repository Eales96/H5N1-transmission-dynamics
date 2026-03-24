# H5N1-transmission-dynamics
Modelling of H5N1 transmission dynamics on dairy farms

## Funding 
This work was funded by the Livestock Biosecurity Fund - Cattle Compensation fund.
This code supports the work done under the agreement CCF 25.19: "Preparing for potential H5N1 incursions into Victorian dairy cattle".

## Reproducing the analyses
This code reproduces the analyses in the technical report produced under milestone 2 of the grant agreement. To reproduce the analyses all 'simulation' scripts should be run, and then all 'figure' scripts:

# Example simulation trajectories
To reproduce the figure showing an example infection trajectory for each transmission regime run the R script 'figures_example_simulation_timeseries.R'.
This will run three simulations (one for each transmission regime) and save the outputs in the 'simulations' subdirectory. It will then create the figure in the technical report. If the simulations have been run previosuly this section of code can be skipped when rerunning.

# Pre-emptive cohorting simulations
To reproduce all analyses of pre-emptive cohorting:
1) Run the R script 'simulations_preemptive.R'. This will save multiple 'parquet' files in the 'simulations' subdirectory.
2) Run the R script 'figures_preemptive.R'. This will reproduce all relevant figures in the technical report.

# Bulk milk sample testing simulations
To reproduce all analyses for bulk milk sample testing:
1) Run the R script 'simulations_testing.R'. This will save multiple 'parquet' files in the 'simulations' subdirectory.
2) Run the R script 'figures_testing.R'. This will reproduce all relevant figures in the technical report.

# Reactive cohorting simulations
To reproduce all analyses of reactive cohorting:
1) Run the R script 'simulations_reactive.R'. This will save multiple 'parquet' files in the 'simulations' subdirectory.
2) Run the R script 'figures_reactive.R'. This will reproduce all relevant figures in the technical report.
