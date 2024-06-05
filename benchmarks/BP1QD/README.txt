The files in this folder can be used to run the community code benchmark problem BP1-QD, defined in Erickson et al. (2020).

The script BP1QD_Plots.m will automatically set up and run the problem in RSFaultZ for a vertical strike-slip fault. BP1QD_Plots.m calls Run_BP1QD.m which sets up the simulation. Then BP1QD_Plots.m runs the simulation and calls BP3QD_OnFault.m, which extracts on-fault simulation data at specified locations.

For comparison with the output from RSFaultZ, this folder also contains data from the BICyclE code, downloaded from  https://strike.scec.org/cvws/cgi-bin/seas.cgi

Erickson et al., (2020). The Community Code Verification Exercise for Simulating Sequences of Earthquakes and Aseismic Slip (SEAS), Seismol. Res. Lett. 91, 874–890, doi: 10.1785/0220190248.
