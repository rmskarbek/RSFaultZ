%%% WARNING: This script will run 70 quasi-static simulations of fault motion, the results
%%% of which are shown shown in Figure 3 in Skarbek, Saffer, and Savage (202X). It would
%%% take a few days of execution time to run all of these simulations on a desktop
%%% computer.


%%%---------------------------------------------------------------------%%%
%%% Simulation parameters and options. See BlocksRun.m for more info.
%%%---------------------------------------------------------------------%%%

%%% Model geometry.
geom = 'Thrust Fault';

%%% Target fault length.
FaultLength = 100e3;

%%% Mean values of (a/b) at 150 C.
aminusb_VS = 0.0054;              % illite/quartz
aminusb_VW = -0.0031;             % carbonate

%%% Options.
calc = 'boundary';
type = 'simulation';
effstressF = 'variable';
symmetry = 'symmetric';

%%%------------------------------------------------------------------------------------%%%
%%% Section 1.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_1 = [0.9; (0.95:0.025:1)'];

%%% Lhat_VW values for simulations.
Lhat_VW_1 = 0.05;

%%% Run the simulations.
[T_1, P_1] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_1, FaultLength, Lhat_VW_1, geom, calc,...
    type, effstressF, symmetry);

%%%------------------------------------------------------------------------------------%%%
%%% Section 2.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_2 = (0.925:0.025:1)';

%%% Lhat_VW values for simulations.
Lhat_VW_2 = (0.1:0.1:0.4)';

%%% Run the simulations.
[T_2, P_2] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_2, FaultLength, Lhat_VW_2, geom, calc,...
    type, effstressF, symmetry);

%%%------------------------------------------------------------------------------------%%%
%%% Section 3.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_3 = 0.95;

%%% Lhat_VW values for simulations.
Lhat_VW_3 = [0.5; 0.55];

%%% Run the simulations.
[T_3, P_3] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_3, FaultLength, Lhat_VW_3, geom, calc,...
    type, effstressF, symmetry);

%%%------------------------------------------------------------------------------------%%%
%%% Section 4.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_4 = 0.975;

%%% Lhat_VW values for simulations.
Lhat_VW_4 = (0.5:0.05:0.65)';

%%% Run the simulations.
[T_4, P_4] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_4, FaultLength, Lhat_VW_4, geom, calc,...
    type, effstressF, symmetry);

%%%------------------------------------------------------------------------------------%%%
%%% Section 5.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_5 = (1:0.05:1.15)';

%%% Lhat_VW values for simulations.
Lhat_VW_5 = (0.5:0.05:0.7)';

%%% Run the simulations.
[T_5, P_5] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_5, FaultLength, Lhat_VW_5, geom, calc,...
    type, effstressF, symmetry);

%%% (a/b)_m values for simulations.
ab_M_5b = 1.05;

%%% Lhat_VW values for simulations.
Lhat_VW_5b = 0.675;

%%% Run the simulations.
[T_5b, P_5b] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_5b, FaultLength, Lhat_VW_5b, geom, calc,...
    type, effstressF, symmetry);

%%%------------------------------------------------------------------------------------%%%
%%% Section 6.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_6 = 1.2;

%%% Lhat_VW values for simulations.
Lhat_VW_6 = (0.55:0.05:0.75)';

%%% Run the simulations.
[T_6, P_6] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_6, FaultLength, Lhat_VW_6, geom, calc,...
    type, effstressF, symmetry);

%%%------------------------------------------------------------------------------------%%%
%%% Section 7.
%%%------------------------------------------------------------------------------------%%%
%%% (a/b)_m values for simulations.
ab_M_7 = (1.25:0.05:1.35)';

%%% Lhat_VW values for simulations.
Lhat_VW_7 = (0.55:0.05:0.8)';

%%% Run the simulations.
[T_7, P_7] = BlocksRun(aminusb_VW, aminusb_VS, ab_M_7, FaultLength, Lhat_VW_7, geom, calc,...
    type, effstressF, symmetry);