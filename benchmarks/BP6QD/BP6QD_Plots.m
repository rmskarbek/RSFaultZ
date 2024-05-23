%%% This script will run the BP6-QD benchmark and then plot the results along with results
%%% from Pierre Dublanchet's code downloaded from: https://strike.scec.org/cvws/cgi-bin/seas.cgi


%%% With the grid spacing set at 25 m, this script executes in ~1 minutes on Rob's computer:
%%% 3.7 GHz processor, 64 GB RAM, Linux OS.

%%% The recommended grid spacing in the banchmark description in 10 m.

%%%--------------------------------------------------------------------------------------%%%
%%%--------------------------------------------------------------------------------------%%%

%%% Set upn the simulation
R = Run_BP6QD;

%%% Run the simulation.
R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
SimData = struct('Params', R.p, 'Output', R.Out);

%%% Close the RSFaultZ instance.
% delete(R);

%%% Create a table of on fault data at specified locations.
OnFaultData = BP6QD_OnFault(SimData);

%%%--------------------------------------------------------------------------------------%%%
%%%--------------------------------------------------------------------------------------%%%
spy = 365.25*24*3600;

%%% On fault data from Dublanchet at 0 km.
C = readtable('Dublanchet_fltst_strk_00.txt');
C5 = readtable('Dublanchet_fltst_strk_05.txt');
% C = readtable('Ozawa_fltst_strk_00b.txt');

%%%--------------------------------------------------------------------------------------%%%
%%% Slip velocity.
figure
subplot(3, 1, 1)
plot(OnFaultData.fltst_strk_p0.Time_s/spy,  OnFaultData.fltst_strk_p0.SlipRate_log10ms, 'k',...
    'LineWidth', 2)
hold on
plot(OnFaultData.fltst_strk_p5.Time_s/spy,  OnFaultData.fltst_strk_p5.SlipRate_log10ms, 'k',...
    'LineWidth', 2)
plot(C.t/spy, C.slip_rate, 'm--', 'LineWidth', 2);
plot(C5.t/spy, C5.slip_rate, 'b--', 'LineWidth', 2);
hold off

xlim([0 2])
ylim([-12 -6])
xlabel('Time (years)')
ylabel('log[Slip Rate (m/s)]')
legend('RSFaultZ', '', 'Dublanchet 0 km', 'Dublanchet +5 km', 'Location', 'northeast');

%%%--------------------------------------------------------------------------------------%%%
%%% Shear stress.
subplot(3, 1, 2)
plot(OnFaultData.fltst_strk_p0.Time_s/spy, OnFaultData.fltst_strk_p0.ShearStress_MPa, 'k',...
    'LineWidth', 2)
hold on
plot(OnFaultData.fltst_strk_p5.Time_s/spy, OnFaultData.fltst_strk_p5.ShearStress_MPa, 'k',...
    'LineWidth', 2)
plot(C.t/spy, C.shear_stress, 'm--', 'LineWidth', 2);
plot(C5.t/spy, C5.shear_stress, 'b--', 'LineWidth', 2);
hold off

xlim([0 2])
ylim([25 31])
xlabel('Time (years)')
ylabel('Shear Stress (MPa)')
legend('RSFaultZ', '', 'Dublanchet 0 km', 'Dublanchet +5 km', 'Location', 'northeast');

%%%--------------------------------------------------------------------------------------%%%
%%% Pore pressure.
subplot(3, 1, 3)
plot(OnFaultData.fltst_strk_p0.Time_s/spy, OnFaultData.fltst_strk_p0.PorePressure_MPa, 'k',...
    'LineWidth', 2)
hold on
plot(OnFaultData.fltst_strk_p5.Time_s/spy, OnFaultData.fltst_strk_p5.PorePressure_MPa, 'k',...
    'LineWidth', 2)
plot(C.t/spy, C.pore_pressure, 'm--', 'LineWidth', 2);
plot(C5.t/spy, C5.pore_pressure, 'm--', 'LineWidth', 2);
hold off

xlim([0 2])
ylim([0 8])
xlabel('Time (years)')
ylabel('Pore Pressure (MPa)')
legend('RSFaultZ', '', 'Dublanchet 0 km', 'Dublanchet +5 km', 'Location', 'northeast');