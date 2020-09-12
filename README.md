# mscprojectminimalwhole-cell
Code from my MSc project, "A Minimal Physiochemical Whole-Cell Model for Engineering Biology; Extending an Existing Bacterial Cell Model". 

Please find below a small description of what each MATLAB script does, and each of these can be found in the folder "Final Code".
Code adapted from Weisse et al., "A mechanistic link between cellular trade-offs, gene expression and growth", PNAS, 2015, and previous masters student's work “The impact of cellular trade-offs on synthetic genetic circuits across cellpopulations", MSc Research Report, 2017.

## Basic Model

cellmodel_odes.m: This file implements the right hand side of the ODE, p and zombie variables removed to match the paper exactly.

cellmodel_Driver.m: This file initializes parameters and calls the solver routine.

cellmodel_odes_PlotsParam.m: This file initializes parameters, calls the solver routine and plots concentrations and transcription rates.

cellmodel_Burden_plot.m: Loops the solver over a range of induction level values, $w_g$, to run repressillator model. 

cellmodel_Monod_Growth.m: Loops the solver over a range of external nutrient values to plot Monod's growth law.


## Population Growth and Dynamic External Nutrient

cellmodel_odes_external.m: This file initializes parameters and calls the solver routine with the inclusion of dynamic external nutrient and population growth. 

cellmodel_Driver_external.m: This file initializes parameters and calls the solver routine with external nutrient and population growth.


## Competing Strains

cellmodel_variednut_2.m: Loops the solver over a range of external nutrient values to plot Monod's growth law, with two strains in a shared environment.


cellmodel_Driver_2.m: This file initializes parameters and calls the solver routine with two strains of cells in the same environment. Used to plot protein concentrations and internal nutrients of both strains.

cellmodel_odes_2.m: This file implements the right hand side of the ODE with external nutrient and population growth with two strains of cells in the same environment.

cellmodel_odes_2_edit.m: This file implements the right hand side of the ODE with two strains of cells in the same environment. Model variables are doubled, except extracellular nutrients to study the dynamics of a competing, more dominant strain.

cellmodel_odes_external_2_edit.m: This file implements the right hand side of the ODE with dynamic external nutrient and population growth with two competitive strains.

cellmodel_PlotParam_external_2.m: This file initializes parameters and calls the solver routine with dynamic external nutrient and population growth, with two competitive strains.
