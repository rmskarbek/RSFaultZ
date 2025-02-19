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

%%%---------------------------------------------------------------------%%%

%%% This script takes ~21 minutes to execute on Rob's computer: 3.7 GHz 
%%% processor, 64 GB RAM, Linux OS.

%%% Edge set with EffStress_0 = 5 MPa took 57 min.

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
ab1 = [0.3; 0.5; 0.7];

%%% Full space. Values of fault length as Length_star = FaultLength/h_star_F.
%%% Note: exp(-1) = 0.3679. 0.36 might be unstable on long enough
%%% timescales; but 0.35 decays, so is definitely stable.
calc = 'length';
Length_hat = [0.35; 0.36; exp(-1); 0.37; 0.3725; 0.375; 0.38; 0.39; 0.4];

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% To re-run the simulations, comment out the load commands on Lines 55
%%% and 56, and uncomment Lines 59 - 61 and 65 - 67.
load('StabilityData_Fig3.mat')
load('StabilityData_Fig3_Edge.mat')

%%% First run with 5 MPa effective stress.
% EffStress_0 = 5;
% [T_L5, P_L5] = SimulationRun1(ab1, Length_hat, geom, beta, mu_0, calc,...
%     Burial_Star, EffStress_0);


%%% Second run with 50 MPa effective stress.
% EffStress_0 = 50;
% [T_L50, P_L50] = SimulationRun1(ab1, Length_hat, geom, beta, mu_0, calc,...
%     Burial_Star, EffStress_0);

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% This loop pulls out the maximum slip speeds for each value of
%%% Length_starF.
V_maxF5 = nan(numel(Length_hat), numel(ab1));
V_maxF50 = nan(numel(Length_hat), numel(ab1));

for i = 1:numel(Length_hat)
    for j = 1:numel(ab1)
        V_maxF5(i, j) = max(T_L5.(i){j,1}.Vmax);
        V_maxF50(i, j) = max(T_L50.(i){j,1}.Vmax);   
    end
end

%%% Get the plate rate, it's the same for each simulation.
v_plate = P_L5.(1){1,1}.Material.PlateRate;

%%%---------------------------------------------------------------------%%%
%%% Review results. Uncomment to view time series plots of the maximum slip velocity
%%% for each simulation.
% spy = 365.25*24*3600;
% for j = 1:numel(ab1)
%     for i = 1:numel(Length_hat)
%         f = figure;
%         semilogy(T_L5.(i){j,1}.Time/spy, T_L5.(i){j,1}.Vmax, 'r')
%         title(sprintf('a/b = %.1f, Lhat = %.4f', [ab1(j), Length_hat(i)]))
%         pause
%         close(f)
%     end
% end
% 
% for j = 1:numel(ab1)
%     for i = 1:numel(Length_hat)
%         f = figure;
%         semilogy(T_L50.(i){j,1}.Time/spy, T_L50.(i){j,1}.Vmax, 'r')
%         title(sprintf('a/b = %.1f, Lhat = %.4f', [ab1(j), Length_hat(i)]))
%         pause
%         close(f)
%     end
% end

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Plot results.
figure
set(gcf,'Color','w')
hold on
box on
cmap = crameri('batlow');
colors = [cmap(256/4,:); cmap(256/2,:); cmap(3*256/4,:)];
P = [{'p_3_5'; 'p_5_5'; 'p_7_5'}, {'p_3_50'; 'p_5_50'; 'p_7_50'}];


for i = 1:numel(ab1)
    eval([P{i,1}, '= plot(Length_hat, V_maxF5(:,i)/v_plate, ''s-'', ''Color'', colors(i,:), ''MarkerFaceColor'', colors(i,:), ''MarkerEdgeColor'', ''k'', ''MarkerSize'', 20, ''LineWidth'', 3);']);
    drawnow
    eval([P{i,1},'.MarkerHandle.LineWidth = 1'])
    % plot(Length_starL, V_maxF5(:,i)/v_plate, 's-', 'Color', colors{i}, 'MarkerFaceColor',...
    %     colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 15)
end

for i = 1:numel(ab1)
    eval([P{i,2}, '= plot(Length_hat, V_maxF50(:,i)/v_plate, ''d-'', ''Color'', colors(i,:), ''MarkerFaceColor'', colors(i,:), ''MarkerEdgeColor'', ''k'', ''MarkerSize'', 15, ''LineWidth'', 3);']);
    drawnow
    eval([P{i,2},'.MarkerHandle.LineWidth = 1'])
    % plot(Length_starL, V_maxF50(:,i)/v_plate, 'd-', 'Color', colors{i}, 'MarkerFaceColor',...
    %     colors{i}, 'MarkerEdgeColor', 'k', 'MarkerSize', 10)
end
plot(exp(-1)*[1 1], [9e-1 2e2], 'k--', 'LineWidth', 3)

%%% Mark locations of the critical fault lengths.
plot(mean(Length_hat(4:5)), 10^mean(log10([V_maxF5(4,1), V_maxF5(5,1)]))/v_plate, 'x',...
    'Color', colors(1,:), 'LineWidth', 3, 'MarkerSize', 15)
plot(mean(Length_hat(5:6)), 10^mean(log10([V_maxF5(5,2), V_maxF5(6,2)]))/v_plate, 'x',...
    'Color', colors(2,:), 'LineWidth', 3, 'MarkerSize', 15)
plot(mean(Length_hat(4:5)), 10^mean(log10([V_maxF5(4,3), V_maxF5(5,3)]))/v_plate, 'x',...
    'Color', colors(3,:), 'LineWidth', 3, 'MarkerSize', 15)

%%% Plot the edge perturbation results.
plot(Length_hat, V_maxF5_Edge(:,1)/v_plate, 'ko', 'MarkerSize', 8)
plot(Length_hat, V_maxF5_Edge(:,2)/v_plate, 'ko', 'MarkerSize', 8)
plot(Length_hat, V_maxF5_Edge(:,3)/v_plate, 'ko', 'MarkerSize', 8)

hold off
xlim([0.349 0.401])
ylim([9e-1 2e2])
ax = gca;
ax.YScale = 'log';
ax.XTick = (0.35:0.01:0.4);
% ax.XTick = [0.35, 0.36, exp(-1), (0.37:0.01:0.4)];
% ax.XTickLabel{3,1} = '$e^{-1}$';
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';

xlabel('Fault Length ($L / h_F^*$)', 'Interpreter', 'latex')
ylabel('Maximum Slip Velocity ($v / v_0$)', 'Interpreter', 'latex')

%%% Add the legend. For tabular legends see:
%%% https://www.mathworks.com/matlabcentral/answers/519942-creating-a-tabular-legend
%%% https://www.mathworks.com/matlabcentral/answers/339686-how-do-i-get-a-tabular-legend
lgd = legend('{ }', '{ }', '{ }', '\ \ \ 0.3', '\ \ \ 0.5', '\ \ \ 0.7', 'Location',...
    'northwest', 'Interpreter', 'latex');
lgd.NumColumns = 2;

title(lgd, ['$L_b$ (km) $\ \ \ \ \ \ $' newline '\ 0.8 \ \ \ \  8 $\ \ \ (a/b)$'],...
    'Interpreter', 'latex')
drawnow

%%% Change the line width of the legend markers.
for i = 1:3
    for j = 1:2
        lineEntry = findobj(lgd.EntryContainer, 'Object', eval(P{i,j}));
        entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
        entryMarker.LineWidth = 1.0;
        % entryMarker.Style = 'none';
    end
end
