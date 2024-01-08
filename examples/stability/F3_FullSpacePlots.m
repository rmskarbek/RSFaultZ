%%%---------------------------------------------------------------------%%%

%%% This script executes one calculation:

%%% Calculation 1: Run full space simulations for different fault lengths.
%%% Use two values of a/b, and two values of L_b. Change L_b by changing
%%% the effective stress to from 5 to 50 MPa. Maximum slip speeds depend on
%%% these values, but the stability boundary does not.

%%%---------------------------------------------------------------------%%%

%%% This script executes the calculations according to the following steps:
%%% 1. Sets up ranges of parameter values in this script.

%%% 2. Calls SimulationRun1.m: executes the simulations for each set of 
%%% parameter values; this occurs in one step.

%%%     2a. SimulationSetup.m: Creates an instance of RSFaultZ; sets the
%%%     parameter values; and runs the simulation.

%%% 3. Generates Figure 3 in Skarbek, 2024 JGR.

%%%---------------------------------------------------------------------%%%

%%% This script takes ~21 minutes to execute on Rob's computer: 3.7 GHz 
%%% processor, 64 GB RAM, Linux OS.

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Calculation 1: Run full space simulations for different fault lengths.
%%% Use two values of a/b, and two values of L_b. Change L_b by changing
%%% the effective stress to from 5 to 50 MPa. Maximum slip speeds depend on
%%% these values, but the stability boundary does not.

%%% Full space geometry does not have a dip angle or burial depth.
geom = 'Full Space';
Burial_Star = nan;
beta = nan;
mu_0 = 0.6;

%%% Values of a/b.
ab1 = [0.7; 0.9];

%%% Full space. Values of fault length as Length_star = FaultLength/h_star_F.
%%% Note: exp(-1) = 0.3679. 0.36 might be unstable on long enough
%%% timescales; but 0.35 decays, so is definitely stable.
calc = 'length';
Length_starL = [0.35; 0.36; exp(-1); 0.37; 0.38; 0.39; (0.4:0.1:1)'];

%%% First run with 5 MPa effective stress.
EffStress_0 = 5;
[T_L5, P_L5] = SimulationRun1(ab1, Length_starL, geom, beta, mu_0, calc,...
    Burial_Star, EffStress_0);

%%% Second run with 50 MPa effective stress.
EffStress_0 = 50;
[T_L50, P_L50] = SimulationRun1(ab1, Length_starL, geom, beta, mu_0, calc,...
    Burial_Star, EffStress_0);

%%% This loop pulls out the maximum slip speeds for each value of
%%% Length_starF.
V_maxF5 = nan(numel(Length_starL), numel(ab1));
V_maxF50 = nan(numel(Length_starL), numel(ab1));

for i = 1:numel(Length_starL)
    for j = 1:numel(ab1)
        V_maxF5(i, j) = max(T_L5.(i){j,1}.Vmax);
        V_maxF50(i, j) = max(T_L50.(i){j,1}.Vmax);
    end
end

%%% Get the plate rate, it's the same for each simulation.
v_plate = P_L5.(1){1,1}.Material.PlateRate;


%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Plot results. MarkerHandle.LineWidth sets the marker border line
%%% thickness separately from the 'LineWidth' setting.
figure;

%%% (a/b) = 0.7, L_b = 0.8 km
p3 = plot(Length_starL, V_maxF50(:,1)/v_plate, 'c--s', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'c', 'MarkerSize', 15, 'LineWidth', 3);
drawnow
p3.MarkerHandle.LineWidth = 1;
hold on;

%%% (a/b) = 0.9, L_b = 0.8 km
p4 = plot(Length_starL, V_maxF50(:,2)/v_plate, 'm--s', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'm', 'MarkerSize', 10, 'LineWidth', 3);
drawnow
p4.MarkerHandle.LineWidth = 1;

%%% (a/b) = 0.7, L_b = 8 km
p1 = plot(Length_starL, V_maxF5(:,1)/v_plate, 'c-s', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'c', 'MarkerSize', 15, 'LineWidth', 3);
drawnow
p1.MarkerHandle.LineWidth = 1;

%%% (a/b) = 0.9, L_b = 8 km
p2 = plot(Length_starL, V_maxF5(:,2)/v_plate, 'm-s', 'MarkerEdgeColor', 'k',...
    'MarkerFaceColor', 'm', 'MarkerSize', 10, 'LineWidth', 3);
drawnow
p2.MarkerHandle.LineWidth = 1;


plot(exp(-1)*[1 1], [1 1e9], 'k--', 'LineWidth', 3)
hold off

xlim([0.3 1])
ylim([1 1e9])

ax = gca;
ax.XTick = [0.3, exp(-1), (0.4:0.1:1)];
ax.XTickLabel{2,1} = '$e^{-1}$';
ax.YScale = 'log';
ax.YTick = logspace(0,9,10);
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';

xlabel('Fault Length ($L / h_F^*$)', 'Interpreter', 'latex')
ylabel('Maximum Slip Velocity ($v / v_0$)', 'Interpreter', 'latex')

%%% Add the legend. For tabular legends see:
%%% https://www.mathworks.com/matlabcentral/answers/519942-creating-a-tabular-legend
%%% https://www.mathworks.com/matlabcentral/answers/339686-how-do-i-get-a-tabular-legend
lgd = legend('{ }', '{ }', '\ \ \ 0.7', '\ \ \ 0.9', 'Location', 'northwest',...
    'Interpreter', 'latex');
lgd.NumColumns = 2;

title(lgd, ['$L_b$ (km) $\ \ \ \ \ \ $' newline '\ 0.8 \ \ \ \  8 $\ \ \ (a/b)$'],...
    'Interpreter', 'latex')
drawnow

%%% Remove the markers from the legend entries.
lineEntry = findobj(lgd.EntryContainer, 'Object',p1);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',p2);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',p3);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',p4);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

% %%%---------------------------------------------------------------------%%%
% %%%---------------------------------------------------------------------%%%
% 
% %%% Calculation 2: Run full space simulations for different perturation 
% %%% lengths. 
% calc = 'perturbation';
% Length_starP = 1./2.^(6:-1:0)';
% [T_P, P_P] = SimulationRun1(ab1(1), Length_starP, geom, beta, mu_0, calc,...
%     Burial_Star);
% 
% %%% This loop pulls out the maximum slip speeds for each value of
% %%% Length_starP.
% V_maxP = nan(numel(Length_starP),1);
% 
% for i = 1:numel(Length_starP)
%     V_maxP(i, 1) = max(T_P.(i){1,1}.Vmax);
% end