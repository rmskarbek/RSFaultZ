%%% This script will run the BP1-QD benchmark and then plot the results 
%%% along with results from BICyclE downloaded from:
%%% https://strike.scec.org/cvws/cgi-bin/seas.cgi

%%% For details, see Erickson et al., (2019), The Community Code 
%%% Verification Exercise for Simulating Sequences of Earthquakes and
%%% Aseismic Slip (SEAS). Seismol. Res. Lett. 91, 874–890, 
%%% doi: 10.1785/0220190248.


%%% This script executes in ~2 minutes on Rob's computer: 3.7 GHz processor,
%%% 64 GB RAM, Linux OS.

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Run 1: Vertical strike slip fault.
DipAngle = 90;

%%% Run the simulation
SimData_SS90 = Run_BP1QD(DipAngle);

%%% Create a table of on fault data at specified locations.
OnFaultData_SS90 = BP3QD_OnFault(SimData_SS90);

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Figures 5C and 5D from Erickson et al. 2019: on-fault shear stress and 
%%% slip rate at 7.5 km along-fault position.
spy = 365.25*24*3600;

%%% On fault data from FDRA at 7.5 km.
C_075_SS90= readtable('Jiang_fltst_dp075_SS90.txt');

figure;
tiledlayout(1,2);

%%% Shear stress.
nexttile;
plot(OnFaultData_SS90.fltst_dp75.Time_s/spy,...
    OnFaultData_SS90.fltst_dp75.ShearStress_MPa, 'k', 'LineWidth', 2)
hold on
plot(C_075_SS90.t/spy, C_075_SS90.shear_stress, 'm--', 'LineWidth', 2);
hold off

xlim([0 700])
ylim([24 48])
xlabel('Time (years)')
ylabel('Shear Stress (MPa)')
title('Vertical Strike-Slip Fault')

l = legend('RSFaultZ', 'BICyclE', 'Location', 'northwest');

%%% Slip rate.
nexttile
plot(OnFaultData_SS90.fltst_dp75.Time_s/spy,...
    OnFaultData_SS90.fltst_dp75.SlipRate_log10ms, 'k', 'LineWidth', 2)
hold on
plot(C_075_SS90.t/spy, C_075_SS90.slip_rate, 'm--', 'LineWidth', 2);
hold off

xlim([0 700])
ylim([-17 1])
xlabel('Time (years)')
ylabel('log[Slip Rate (m/s)]')
title('Vertical Strike-Slip Fault')