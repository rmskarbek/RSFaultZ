The files in this folder can be used to run the community code benchmark problem BP3-QD, described at https://strike.scec.org/cvws/seas/benchmark_descriptions.html

The script BP6QD_Plots.m will automatically set up and run the problem in RSFaultZ. BP6QD_Plots.m calls Run_BP6QD.m which sets up the simulation. Then BP6QD_Plots.m runs the simulation and calls BP6QD_OnFault.m, which extracts on-fault simulation data at specified locations.
