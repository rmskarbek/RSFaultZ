%%% This script conducts a numerical linear stability analysis for a dipping thrust fault
%%% with heterogeneous friction properties. These results are displayed in Figure 3 in
%%% Skarbek, Saffer, and Savage (202X).

%%% This script executes in about 20 minutes on Rob's computer:  3.7 GHz processor, 64 GB
%%% RAM, Linux OS.

%%%------------------------------------------------------------------------------------%%%
%%% Parameters and options for linear stability calculation. See BlocksRun.m for more info.
%%%------------------------------------------------------------------------------------%%%

%%% Model geometry.
geom = 'Thrust Fault';

%%% Target fault length.
FaultLength = 100e3;

%%% Mean values of (a/b) at 150 C.
aminusb_VS = 0.0054;              % illite/quartz
aminusb_VW = -0.0031;             % carbonate

%%% Options.
calc = 'boundary';
type = 'stability';
effstressF = 'variable';
symmetry = 'symmetric';

%%% Values of (a/b)_m at which to find the partial stability boundary.
AB_M = (1.0125:0.0125:1.3625)';

%%% A dummy value of Lhat_VW is needed as input, but this value is not used for the
%%% stability calculation
Lhat_VW = 0;

%%% Carry out the stability calculation;
[T, P] = BlocksRun(aminusb_VW, aminusb_VS, AB_M, FaultLength, Lhat_VW, geom, calc, type,...
    effstressF, symmetry);

%%% Get the critical VW block sizes.
L_star_VW = cell2mat(T.Lhat_00);

%%% Get the constructed values of (a/b)_m.
AB_m_p = nan(size(L_star_VW));
for i = 1:numel(P)
    AB_m_p(i,1) = P{i,1}{1,1}.Friction.ab_M;
end

%%% Save the results.
save('PartialStabilityData.mat')