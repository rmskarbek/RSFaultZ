The files in this folder can be used to run the community code benchmark problem BP1-QD, defined in Erickson et al. (2019).

The script BP1QD_Plots.m will automatically set up and run the problem in RSFaultZ for a vertical strike-slip fault. BP1QD_Plots.m calls Run_BP1QD.m which executes the simulation, then calls BP3QD_OnFault.m, which extracts on-fault simulation data at specified locations.

For comparison with the output from RSFaultZ, this folder also contains data from the BICyclE code, downloaded from  https://strike.scec.org/cvws/cgi-bin/seas.cgi

Erickson et al., (2023), Incorporating Full Elastodynamic Effects and Dipping Fault Geometries in Community Code Verification Exercises for Simulations of Earthquake Sequences and  Aseismic Slip (SEAS), Bull. Seismol. Soc. Am. 113, 499–523, doi: 10.1785/0120220066.
