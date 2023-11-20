# CenoCO2
Code and data for Baysian timeseries models of Cenozoic pCO<sub>2</sub> and temperature published by the CenCO<sub>2</sub>PIP Consortium.

## data/
Data used in analyses.
- **230902_proxies.xlsx** Multiproxy paleo-CO<sub>2</sub> data
- **Hansen13Russell.csv** CO<sub>2</sub> reconstruction from Hansen, et al. (2013)
- **Westerhold.xlsx** CENOGRID benthic &delta;<sup>18</sup>O data and paleotemperatures from Westerhold, et al. (2020)

## code/
Scripts used for data analysis and plotting. Some output is saved to folder bigout/ which is not included in archive; this output can be recreated by running the scripts.
- **1_Driver.R** Load data and run CO<sub>2</sub> inversion
- **2_DriverT.R** Load data and run temperature inversion
- **3_DignosticsStats.R** Review inversion metrics and output, also prep paleotemperature dataset from Ring, et al. (2023)
- **4_MainPlots.R** Plotting for main manuscript
- **5_SIPlots.R** Supplemental information plotting
- **6_OtherPlots.R** Plotting for presentations and press
- **Helpers.R** Functions used in other scripts
- **PrepForPlots.R** Functions used in other scripts

## code/models/
JAGS models used to run inversions
- **model.R** Model for paleo-CO<sub>2</sub>
- **model_T.R** Model for paleotemperature

## out/
Output generated by scripts.
- **1MyrCO2.csv** Quantile values for ln(CO<sub>2</sub>) reconstructed at 1 Myr resolution
- **110kryCO2.csv** Quantile values for ln(CO<sub>2</sub>) reconstructed at 100 kyr resolution
- **500kyrCO2.csv** Quantile values for ln(CO<sub>2</sub>) reconstructed at 500 kyr resolution
- **500kyrCO2MarOnly.csv** Quantile values for ln(CO<sub>2</sub>) reconstructed at 500 kyr resolution using only marine proxies
- **500kyrTemp.csv** Quantile values for temperature reconstructed at 500 kyr resolution
- **RingTemp.csv** Quantile values for Ring, et al. (2023) temperature reconstruction and associated ln(CO<sub>2</sub>)

## out/main_figs/
Figures and figure components included in the main manuscript.
- **Fig2.eps** Components used in Figure 2
- **Fig3.png** Figure 3
- **PrintFig.pdf** Warming stripes figure included in print summary

## out/SI_figs/
Figures included as supplementary information. File names correspond to figure numbers in the supplementary text.

## out/other_figs/
Figures generated for presentations and press.
- **Age_uncert.png** Figure illustrating the incorporation and impact of age model uncertainty in the paleo-CO<sub>2</sub> inversion
- **CenozoicCO2_slide.png** Simple depiction of 500 kyr resoltution CO<sub>2</sub> and temperature reconstructions
- **DataDens1.png** Shows the number of proxy data points represented within 1 and 5 Myr moving windows across the Cenozoic
- **DataDens2.png** Shows the number of unique proxy types represented within 1 and 5 Myr moving windows across the Cenozoic
- **Lamont_press.png** Warming stripes figure version prepared for Lamont press office
- **UU_press.png** Warming stripes figure version prepared for University of Utah press office