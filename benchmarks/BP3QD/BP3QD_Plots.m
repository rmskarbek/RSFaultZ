%%% This script will run two versions of the BP3-QD benchmark and then plot
%%% the results along with results for the FDRA model downloaded from:
%%% https://strike.scec.org/cvws/cgi-bin/seas.cgi

%%% For details, see Erickson et al., (2023), Incorporating Full 
%%% Elastodynamic Effects and Dipping Fault Geometries in Community Code
%%% Verification Exercises for Simulations of Earthquake Sequences and 
%%% Aseismic Slip (SEAS), Bull. Seismol. Soc. Am. 113, 499–523, 
%%% doi: 10.1785/0120220066.

%%% This script executes in ~8 minutes on Rob's computer: 3.7 GHz processor,
%%% 64 GB RAM, Linux OS.

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Run 1: Thrust fault with 60 deg dip.
geom = 'Thrust Fault';
DipAngle = 60;

%%% Run the simulation
SimData_T60 = Run_BP3QD(geom, DipAngle);

%%% Create a table of on fault data at specified locations.
OnFaultData_T60 = BP3QD_OnFault(SimData_T60);

%%% Create tables of displacement on the upper surface at specified
%%% locations.
[C, Time_T60, T_UX_T60, T_UY_T60] = BP3QD_Displacements(SimData_T60);


%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Run 2: Normal fault with 30 deg dip.
geom = 'Normal Fault';
DipAngle = 30;

%%% Run the simulation
SimData_N30 = Run_BP3QD(geom, DipAngle);

%%% Create a table of on fault data at specified locations.
OnFaultData_N30 = BP3QD_OnFault(SimData_N30);

%%% Create tables of displacement on the upper surface at specified
%%% locations.
[C, Time_N30, T_UX_N30, T_UY_N30] = BP3QD_Displacements(SimData_N30);


%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Figure 9 from Erickson et al. 2023: on-fault shear stress and slip
%%% rate at 7.5 km along-fault position.
spy = 365.25*24*3600;

%%% On fault data from FDRA at 7.5 km for both thrust and normal faults.
%%% Thrust fault.
C_075_T60= readtable('cattania_fltst_dp075_T60.txt');
%%% Normal fault.
C_075_N30= readtable('cattania_fltst_dp075_N30.txt');

figure;
tiledlayout(2,2);

%%% Shear stress for thrust fault.
nexttile;
plot(OnFaultData_T60.fltst_dp75.Time_s/spy,...
    OnFaultData_T60.fltst_dp75.ShearStress_MPa, 'k', 'LineWidth', 2)
hold on
plot(C_075_T60.t/spy, C_075_T60.shear_stress, 'm--', 'LineWidth', 2);
hold off

xlim([0 700])
ylim([25 50])
ylabel('Shear Stress (MPa)')
title('60\circ thrust')

l = legend('RSFaultZ', 'FDRA', 'Location', 'northwest');

%%% Slip rate for thrust fault.
nexttile
plot(OnFaultData_T60.fltst_dp75.Time_s/spy,...
    OnFaultData_T60.fltst_dp75.SlipRate_log10ms, 'k', 'LineWidth', 2)
hold on
plot(C_075_T60.t/spy, C_075_T60.slip_rate, 'm--', 'LineWidth', 2);
hold off

xlim([0 700])
ylim([-18 2])
ylabel('log[Slip Rate (m/s)]')
title('60\circ thrust')

%%% Shear stress for normal fault.
nexttile;
plot(OnFaultData_N30.fltst_dp75.Time_s/spy,...
    -OnFaultData_N30.fltst_dp75.ShearStress_MPa, 'k', 'LineWidth', 2)
hold on
plot(C_075_N30.t/spy, C_075_N30.shear_stress, 'm--', 'LineWidth', 2);
hold off

xlim([0 700])
ylim([-50 -25])
xlabel('Time (years)')
ylabel('Shear Stress (MPa)')
title('30\circ normal')

%%% Slip rate for normal fault.
nexttile
plot(OnFaultData_N30.fltst_dp75.Time_s/spy,...
    OnFaultData_N30.fltst_dp75.SlipRate_log10ms, 'k', 'LineWidth', 2)
hold on
plot(C_075_N30.t/spy, C_075_N30.slip_rate, 'm--', 'LineWidth', 2);
hold off

xlim([0 700])
ylim([-18 2])
xlabel('Time (years)')
ylabel('log[Slip Rate (m/s)]')
title('30\circ normal')


%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Figure 12 from Erickson et al., 2023; but just using the 60 deg thrust
%%% and 30 deg normal fault simulations here. Vertical and horizontal
%%% surface displacements at p\m 16 km and p\m 0 km.

%%% Surface data from FDRA for +16 km.
C_p16_T60 = readtable('cattania_srfst_fn_p16_T60.txt');
C_p16_N30 = readtable('cattania_srfst_fn_p16_N30.txt');

%%% Surface data for FDRA -16 km.
C_n16_T60 = readtable('cattania_srfst_fn_n16_T60.txt');
C_n16_N30 = readtable('cattania_srfst_fn_n16_N30.txt');

%%% Surface data for FDRA +0 km.
C_p0_T60 = readtable('cattania_srfst_fn_p0_T60.txt');
C_p0_N30 = readtable('cattania_srfst_fn_p0_N30.txt');

%%% Surface data from FDRA for -0 km is not provided at:
%%% https://strike.scec.org/cvws/cgi-bin/seas.cgi

figure;
tiledlayout(2,2);

%%% Thrust fault horizontal displacements.
nexttile
plot(Time_T60/spy, T_UX_T60.p_16km, 'k', 'LineWidth', 2)
hold on
plot(Time_T60/spy, T_UX_T60{:,5}, 'k', 'LineWidth', 2)
plot(Time_T60/spy, T_UX_T60{:,6}, 'k', 'LineWidth', 2)
plot(Time_T60/spy, T_UX_T60.n_16km, 'k', 'LineWidth', 2)

plot(C_p16_T60.t/spy, C_p16_T60.disp_1, 'm--', 'LineWidth', 2)
plot(C_n16_T60.t/spy, C_n16_T60.disp_1, 'm--', 'LineWidth', 2)
plot(C_p0_T60.t/spy, C_p0_T60.disp_1, 'm--', 'LineWidth', 2)

hold off

xlim([0 600])
ylim([-10 10])
ylabel('Horizontal Displacement (m)')
title('60\circ thrust')

ax = gca;
ax.YTick = -10:2:10;

l = legend('RSFaultZ', '', '', '', 'FDRA', 'Location', 'northwest');

%%% Normal fault horizontal displacements.
nexttile
plot(Time_N30/spy, T_UX_N30.p_16km, 'k', 'LineWidth', 2)
hold on
plot(Time_N30/spy, T_UX_N30{:,5}, 'k', 'LineWidth', 2)
plot(Time_N30/spy, T_UX_N30{:,6}, 'k', 'LineWidth', 2)
plot(Time_N30/spy, T_UX_N30.n_16km, 'k', 'LineWidth', 2)

plot(C_p16_N30.t/spy, C_p16_N30.disp_1, 'm--', 'LineWidth', 2)
plot(C_n16_N30.t/spy, C_n16_N30.disp_1, 'm--', 'LineWidth', 2)
plot(C_p0_N30.t/spy, C_p0_N30.disp_1, 'm--', 'LineWidth', 2)
hold off

xlim([0 600])
ylim([-10 10])
title('30\circ normal')

ax = gca;
ax.YTick = -10:2:10;

%%% Thrust fault vertical displacements.
nexttile
plot(Time_T60/spy, T_UY_T60.p_16km, 'k', 'LineWidth', 2)
hold on
plot(Time_T60/spy, T_UY_T60{:,5}, 'k', 'LineWidth', 2)
plot(Time_T60/spy, T_UY_T60{:,6}, 'k', 'LineWidth', 2)
plot(Time_T60/spy, T_UY_T60.n_16km, 'k', 'LineWidth', 2)

plot(C_p16_T60.t/spy, C_p16_T60.disp_2, 'm--','LineWidth', 2)
plot(C_n16_T60.t/spy, C_n16_T60.disp_2, 'm--','LineWidth', 2)
plot(C_p0_T60.t/spy, C_p0_T60.disp_2, 'm--','LineWidth', 2)
hold off

xlim([0 600])
ylim([-10 10])
xlabel('Time (years)')
ylabel('Vertical Displacement (m)')
title('60\circ thrust')

ax = gca;
ax.YTick = -10:2:10;

%%% Normal fault vertical displacements.
nexttile
plot(Time_N30/spy, T_UY_N30.p_16km, 'k', 'LineWidth', 2)
hold on
plot(Time_N30/spy, T_UY_N30{:,5}, 'k', 'LineWidth', 2)
plot(Time_N30/spy, T_UY_N30{:,6}, 'k', 'LineWidth', 2)
plot(Time_N30/spy, T_UY_N30.n_16km, 'k', 'LineWidth', 2)

plot(C_p16_N30.t/spy, C_p16_N30.disp_2, 'm--', 'LineWidth', 2)
plot(C_n16_N30.t/spy, C_n16_N30.disp_2, 'm--', 'LineWidth', 2)
plot(C_p0_N30.t/spy, C_p0_N30.disp_2, 'm--', 'LineWidth', 2)
hold off

xlim([0 600])
ylim([-10 10])
xlabel('Time (years)')
title('30\circ normal')

ax = gca;
ax.YTick = -10:2:10;