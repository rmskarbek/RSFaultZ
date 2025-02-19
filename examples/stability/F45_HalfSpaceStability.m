%%%---------------------------------------------------------------------%%%

%%% This script executes three calculations:

%%% Calculation 1 and 2: Numerically determine the critical fault length 
%%% for dipping faults, by finding the smallest fault length that is 
%%% unstable, for different values of dip angle, burial depth, and friction 
%%% coefficient. This is done for both thrust and normal faults.

%%% This script executes the calculations according to the following steps:
%%% 1. Sets up ranges of parameter values in this script.

%%% 2. Calls StabilityRun1.m: executes the stability analysis for each set
%%% of parameter values; this occurs in two steps.

%%%     2a. StabilitySetup.m: Creates an instance of RSFaultZ and sets the
%%%     parameter values.

%%%     2b. StabilityAnalysis.m: Executes the numerical stability analysis
%%%     using the parameters structure that is created by StabilitySetup.m


%%%---------------------------------------------------------------------%%%

%%% Calculation 3: Computes normalized stressing rates on a dipping fault 
%%% with a uniform velocity distribution. This is carried out by calling
%%% StressingRates.m

%%%---------------------------------------------------------------------%%%

%%% This script takes ~3.3 hours to execute on Rob's computer: 3.7 GHz 
%%% processor, 64 GB RAM, Linux OS.

%%% In total, 4700 individual numerical stability calculations are executed.
%%% The long computation time is due to this, and the small grid spacing:
%%% dxi = H_s_F/250 or dxi = L_b/80, whichever is smaller.


%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Plots the critical fault lengths as a function of dip angle.

%%% Set the geometry for these calculations.
geom = [{'Thrust Fault'}; {'Normal Fault'}];
calc = 'boundary';
type = 'search';

%%% Values of fault dip and normalized burial depth for the calculations.
beta = (0:0.5:90)';
% beta = [0 0.5 10 90]';
Burial_hat = [0; 0.001; 0.01; 0.1; 1];

%%% Values of a/b.
ab = 0.8;

%%% StabilityRun1 does not use the value of Length_hat for 'Thrust Fault'
%%% and 'Normal Fault' geometries.
Length_hat = nan; % dummy value.

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% To re-run the simulations, comment out the load command on Lines 65
%%% and uncomment Lines 68 - 93.
load('StabilityData_Fig4.mat')

%%% Thrust faults.
% j = 1;
% mu_0 = 0.6;
% [T_T06, P_T06] = StabilityRun1(ab, Length_hat, calc, geom{j}, beta, mu_0,...
%     Burial_hat, type);
% 
% mu_0 = 0.2;
% [T_T02, P_T02] = StabilityRun1(ab, Length_hat, calc, geom{j}, beta, mu_0,...
%     Burial_hat(1:4), type);
% 
% mu_0 = 1;
% [T_T1, P_T1] = StabilityRun1(ab, Length_hat, calc, geom{j}, beta, mu_0,...
%     Burial_hat(1:4), type);
% 
% %%% Normal faults.
% j = 2;
% mu_0 = 0.6;
% [T_N06, P_N06] = StabilityRun1(ab, Length_hat, calc, geom{j}, beta, mu_0,...
%     Burial_hat, type);
% 
% mu_0 = 0.2;
% [T_N02, P_N02] = StabilityRun1(ab, Length_hat, calc, geom{j}, beta, mu_0,...
%     Burial_hat(1:4), type);
% 
% mu_0 = 1;
% [T_N1, P_N1] = StabilityRun1(ab, Length_hat, calc, geom{j}, beta, mu_0,...
%     Burial_hat(1:4), type);

%%%-------------------------------------------------------------------------------%%%
%%%-------------------------------------------------------------------------------%%%
figure;
set(gcf,'Color','w')
layout1 = tiledlayout(1, 2, 'Padding', 'tight');
col = {'m'; 'b'; 'c'; 'g'; 'k'};
nexttile(layout1)
hold on;
box on;

%%% Plot thrust fault results.
for i = 1:numel(Burial_hat)
    if i == 1
        plot(beta(2:end), cell2mat(T_T06.(i)), col{i}, 'LineWidth', 3);
    else
        plot(beta, cell2mat(T_T06.(i)), col{i}, 'LineWidth', 3);
    end
end

for i = 1:numel(Burial_hat(1:4))
    if i == 1
        plot(beta(2:end), cell2mat(T_T02.(i)), col{i}, 'LineWidth', 3,...
            'LineStyle', '--');
        plot(beta(2:end), cell2mat(T_T06.(i)), col{i}, 'LineWidth', 3);
        plot(beta(2:end), cell2mat(T_T1.(i)), col{i}, 'LineWidth', 3,...
            'LineStyle', ':');
    else
        plot(beta, cell2mat(T_T02.(i)), col{i}, 'LineWidth', 3,...
            'LineStyle', '--');
        plot(beta, cell2mat(T_T1.(i)), col{i}, 'LineWidth', 3,...
            'LineStyle', ':');
    end
end
hold off

%%% Labels and such.
title('Thrust Faults', 'Interpreter', 'latex', 'FontSize', 20)
text(0.05, 0.05, 'A', 'units', 'normalized', 'FontSize', 30, 'Interpreter',...
    'latex')
xlabel('Dip Angle $\beta$ (deg)', 'Interpreter', 'latex')
ylabel('Critical Fault Length ($L^*_D / h_F^*$)', 'Interpreter', 'latex') 
xlim([0 90])
ylim([0 0.4])
ax = gca;
ax.FontSize = 22;
ax.XTick = 0:10:90;
ax.YTick = [(0:0.05:0.35) exp(-1) 0.4];
ax.YTickLabel{9,1} = '$e^{-1}$';
ax.TickLabelInterpreter = 'latex';

%%%---------------------------------------------------------------------%%%
%%% Plot normal fault results.
nexttile(layout1)
hold on;
box on;
for i = 1:numel(Burial_hat)
    if i == 1
        plot(beta(2:end), cell2mat(T_N06.(i)), col{i}, 'LineWidth', 3);
    else
        plot(beta, cell2mat(T_N06.(i)), col{i}, 'LineWidth', 3);
    end
end

for i = 1:numel(Burial_hat(1:4))
    if i == 1
        m02 = plot(beta(2:end), cell2mat(T_N02.(i)), col{i}, 'LineWidth',...
            3, 'LineStyle', '--');
        m06 = plot(beta(2:end), cell2mat(T_N06.(i)), col{i}, 'LineWidth',...
            3);
        m1 = plot(beta(2:end), cell2mat(T_N1.(i)), col{i}, 'LineWidth',...
            3, 'LineStyle', ':');
    else
        plot(beta, cell2mat(T_N02.(i)), col{i}, 'LineWidth', 3,...
            'LineStyle', '--');
        plot(beta, cell2mat(T_N1.(i)), col{i}, 'LineWidth', 3,...
            'LineStyle', ':');
    end
end
hold off

%%% Labels and such.
title('Normal Faults', 'Interpreter', 'latex', 'FontSize', 20)
text(0.05, 0.05, 'B', 'units', 'normalized', 'FontSize', 30, 'Interpreter',...
    'latex')
xlabel('Dip Angle $\beta$ (deg)', 'Interpreter', 'latex')
ylabel('Critical Fault Length ($L^*_D / h_F^*$)', 'Interpreter', 'latex') 
xlim([0 90])
ylim([0 0.4])
ax = gca;
ax.FontSize = 22;
ax.XTick = 0:10:90;
ax.YTick = [(0:0.05:0.35) exp(-1) 0.4];
ax.YTickLabel{9,1} = '$e^{-1}$';
ax.TickLabelInterpreter = 'latex';

lgd = legend('0', '0.001', '0.01', '0.1', '1', 'Location', 'southeast',...
    'Interpreter', 'latex');
title(lgd, '$d / h^*_F$', 'Interpreter', 'latex')

ax2 = axes('position', ax.Position, 'visible', 'off');
lgd2 = legend(ax2, [m02 m06 m1], '0.2', '0.6', '1', 'Location', 'south',...
    'Interpreter', 'latex');
lgd2.FontSize = 18;
%lgd2.Position = [0.7292    0.1800    0.0740    0.1770];
title(lgd2, '$\mu_0$', 'Interpreter', 'latex')

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Plot critical fault lengths for a few values of dip angle, as a
%%% function of burial depth.
[~, i_05] = min(abs(beta - 0.5));
[~, i_10] = min(abs(beta - 10));
[~, i_20] = min(abs(beta - 20));

figure;
set(gcf,'Color','w')
layout2 = tiledlayout(1, 2, 'Padding', 'tight');
nexttile(layout2)

%%% Thrust faults.
ct = [217 83 25]/255;
t1 = plot(Burial_hat, cell2mat(T_T06{i_05,:})', ':o', 'Color', ct,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', ct, 'MarkerSize', 10,...
    'LineWidth', 3);
% drawnow
% t1.MarkerHandle.LineWidth = 1;
hold on
t2 = plot(Burial_hat, cell2mat(T_T06{i_10,:})', '-o', 'Color', ct,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', ct, 'MarkerSize', 10,...
    'LineWidth', 3);
% drawnow
% t2.MarkerHandle.LineWidth = 1;
t3 = plot(Burial_hat, cell2mat(T_T06{i_20,:})', '--o', 'Color', ct,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', ct, 'MarkerSize', 10,...
    'LineWidth', 3);
% drawnow
% t3.MarkerHandle.LineWidth = 1;


%%% Normal faults.
cn = [114 189 255]/255;
n1 = plot(Burial_hat, cell2mat(T_N06{i_05,:})', ':s', 'Color', cn,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', cn, 'MarkerSize', 10,...
    'LineWidth', 3);
% drawnow
% n1.MarkerHandle.LineWidth = 1;
n2 = plot(Burial_hat, cell2mat(T_N06{i_10,:})', '-s', 'Color', cn,...
    'MarkerEdgeColor', 'k',    'MarkerFaceColor', cn, 'MarkerSize', 10,...
    'LineWidth', 3);
% drawnow
% n2.MarkerHandle.LineWidth = 1;
n3 = plot(Burial_hat, cell2mat(T_N06{i_20,:})', '--s', 'Color', cn,...
    'MarkerEdgeColor', 'k', 'MarkerFaceColor', cn, 'MarkerSize', 10,...
    'LineWidth', 3);
% drawnow
% n3.MarkerHandle.LineWidth = 1;

%%% exp(-1) line.
plot([1e-3 1], exp(-1)*[1 1], 'k', 'LineWidth', 3)
hold off

xlim([1e-3 1])
ylim([0 0.4])
ax = gca;
ax.FontSize = 20;
ax.XScale = 'log';
ax.XTick = [1e-3 1e-2 1e-1 1];
ax.YTick = [(0:0.05:0.35) exp(-1) 0.4];
ax.YTickLabel{9,1} = '$e^{-1}$';
ax.TickLabelInterpreter = 'latex';

xlabel('Normalized Burial Depth ($d / h^*_F$)', 'Interpreter', 'latex')
ylabel('Critical Fault Length ($L^*_D / h_F^*$)', 'Interpreter', 'latex')

text(0.05, 0.05, 'A', 'units', 'normalized', 'FontSize', 30,...
    'Interpreter', 'latex')

lgd = legend('0.5', '10', '20', '0.5', '10', '20', 'Location', 'southeast',...
'Interpreter', 'latex');
lgd.NumColumns = 2;
title(lgd, '$\beta$ (deg)', 'Interpreter', 'latex')

drawnow
t1.MarkerHandle.LineWidth = 1;
t2.MarkerHandle.LineWidth = 1;
t3.MarkerHandle.LineWidth = 1;
n1.MarkerHandle.LineWidth = 1;
n2.MarkerHandle.LineWidth = 1;
n3.MarkerHandle.LineWidth = 1;

% ax2 = axes('position', ax.Position, 'visible', 'off');
% lgd2 = legend(ax2, [t2b n2], 'thrust', 'normal', '1', 'Location', 'south',...
%     'Interpreter', 'latex');
% lgd2.FontSize = 18;

% drawnow

lineEntry = findobj(lgd.EntryContainer, 'Object',t1);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',t2);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',t3);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',n1);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',n2);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

lineEntry = findobj(lgd.EntryContainer, 'Object',n3);
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.Style = 'none';

% lineEntry = findobj(lgd2.EntryContainer, 'Object',n2);
% entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
% entryMarker.Style = 'none';

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Compute normalized stressing rates for a uniform velocity distribution.
[Beta, Shear, G_thrust, G_normal] = StressingRates;

%%% Plot the results.
nexttile(layout2)
plot(Beta, abs(Shear), 'k', 'LineWidth', 3);
hold on;

ct = [217 83 25]/255;
plot(Beta, G_thrust(:,1), '--', 'Color', ct, 'LineWidth', 3);
plot(Beta, G_thrust(:,2), '-', 'Color', ct, 'LineWidth', 3);
plot(Beta, G_thrust(:,3), ':', 'Color', ct, 'LineWidth', 3);

cn = [114 189 255]/255;
plot(Beta, G_normal(:,1), '--', 'Color', cn, 'LineWidth', 3);
plot(Beta, G_normal(:,2), '-', 'Color', cn, 'LineWidth', 3);
plot(Beta, G_normal(:,3), ':', 'Color', cn, 'LineWidth', 3);

tB1 = plot(Beta, G_thrust(:,2), '-', 'Color', ct, 'LineWidth', 3);
nB1 = plot(Beta, G_normal(:,2), '-', 'Color', cn, 'LineWidth', 3);
hold off

xlabel('Dip Angle $\beta$ (deg)', 'Interpreter', 'latex')
ytext = ['Normalized Stressing Rate' newline ...
    '$(2 \pi L)(\dot{\tau} - \mu_0 \dot{\sigma})/(G'' \Delta v)$'];
ylabel(ytext, 'Interpreter', 'latex');
xlim([0 90])
ylim([0 3.5])
ax = gca;
ax.FontSize = 20;
ax.XTick = 0:10:90;    
ax.TickLabelInterpreter = 'latex';

text(0.05, 0.05, 'B', 'units', 'normalized', 'FontSize', 30,...
    'Interpreter', 'latex')

lgdB = legend('0', '0.2', '0.6', '1', '0.2', '0.6', '1', 'Location', 'southeast',...
'Interpreter', 'latex');
title(lgdB, '$\mu_0$', 'Interpreter', 'latex')
lgdB.NumColumns = 2;

axB2 = axes('position', ax.Position, 'visible', 'off');
lgdB2 = legend(axB2, [tB1 nB1], 'thrust', 'normal', '1', 'Location', 'south',...
    'Interpreter', 'latex');
lgdB2.FontSize = 18;