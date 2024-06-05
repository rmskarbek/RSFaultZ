The files in this folder can be used to run the community code benchmark problem BP3-QD, defined in Erickson et al. (2023).

The script BP3QD_Plots.m will automatically set up and run the problem in RSFaultZ for: (1) a thrust fault with 60 deg dip, and (2) normal fault with 30 deg dip. BP3QD_Plots.m calls Run_BP3QD.m which sets up the simulations. Then BP3QD_Plots.m runs the simulation and calls BP3QD_OnFault.m, which extracts on-fault simulation data at specified locations. Finally, BP3QD_Plots.m calls BP3QD_Displacements.m, which computes displacemts at specified locations off of the fault.

For comparison with the output from RSFaultZ, this folder also contains data from the FDRA code, downloaded from  https://strike.scec.org/cvws/cgi-bin/seas.cgi

Erickson et al., (2023), Incorporating Full Elastodynamic Effects and Dipping Fault Geometries in Community Code Verification Exercises for Simulations of Earthquake Sequences and  Aseismic Slip (SEAS), Bull. Seismol. Soc. Am. 113, 499–523, doi: 10.1785/0120220066.
