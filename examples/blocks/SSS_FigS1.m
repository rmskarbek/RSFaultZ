%%% This script generates Figure S1 from Skarbek, Saffer, and Savage (202X). The figure
%%% shows how numerically constructed fault lengths depend on the values of (a/b)_m and
%%% the velocity-weakening block length.

%%% This script generates the figure by generating the numerical grids, rather than by
%%% loading previously generated data.

%%% This script should execute in less than one second.

%%%-------------------------------------------------------------------------------%%%
%%% Carry out the grid calculation to get the constructed fault lengths.
%%%-------------------------------------------------------------------------------%%%

%%% Model geometry.
geom = 'Thrust Fault';

%%% Target fault length.
FaultLength = 100e3;

%%% Mean values of (a/b) at 150 C, needed for grid construction.
aminusb_VS = 0.0054;              % illite/quartz
aminusb_VW = -0.0031;             % carbonate

%%% Options.
calc = 'boundary';
effstressF = 'variable';
symmetry = 'symmetric';

%%% (a/b)_m values for grids.
ab_M = (0.8:0.05:1.5)';

%%% Normalized velocity-weakening block length values for grids.
Lhat_VW = (0.05:0.05:1)';

%%% Compute the grids.
type = 'grid';
[~, G] = BlocksRun(aminusb_VW, aminusb_VS, ab_M, FaultLength, Lhat_VW, geom, calc,...
    type, effstressF, symmetry);

%%% Results.
FL = nan(numel(ab_M), numel(Lhat_VW));
LH = nan(numel(ab_M), numel(Lhat_VW));
AB_M = nan(numel(ab_M), numel(Lhat_VW));
Eta = nan(numel(ab_M), numel(Lhat_VW));
for i = 1:numel(ab_M)           % rows
    for j = 1:numel(Lhat_VW)    % columns
        FL(i, j) = G.(j){i,1}.FaultLength;
        LH(i, j) = G.(j){i,1}.Lhat_VW;
        AB_M(i, j) = G.(j){i,1}.ab_M;
        Eta(i, j) = G.(j){i,1}.eta;
    end
end
%%% Reshape into vectors for scatter plots.
N = numel(FL);
FL = reshape(FL, [N 1]);
LH = reshape(LH, [N 1]);
AB_M = reshape(AB_M, [N 1]);
Eta = reshape(Eta, [N 1]);

%%%-------------------------------------------------------------------------------%%%
%%% Create a scatter plot of the constructed fault lengths.
%%%-------------------------------------------------------------------------------%%%
%%% Set up the figure.
font = 'Palatino Linotype';
figure('DefaultTextFontName',font,'DefaultAxesFontName',font, 'DefaultAxesFontSize', 16,...
    'Units', 'inches', 'Position', [1 1 8.5 6])
set(gcf,'Color','w')

%%% Plot the fault lengths;
scatter(LH, AB_M, 200, FL/100e3, 'filled', 'MarkerEdgeColor', 'k')
box on
grid on
hold on

%%% Use different color schemes for fault lengths outside of +/- 10% range.
cmap1 = crameri('buda');
cmap2 = crameri('batlow');
cmap3 = crameri('hawaii');

i_map90 = round(((90e3 - min(FL))/(max(FL) - min(FL)))*size(cmap2,1));
cmap1(1:i_map90+1,:) = cmap2(1:i_map90+1,:);

i_map110 = round(((max(FL) - 110e3)/(max(FL) - min(FL)))*size(cmap2,1));
i_map110 = size(cmap2,1) - i_map110;
cmap1(i_map110+1:end,:) = cmap3(i_map110+1:end,:);

colormap(cmap1);
C = colorbar('TickLabelInterpreter', 'latex');
C.Location = 'northoutside';
C.Title.String = 'Normalized Constructed Fault Length ($L_c / L$)';
C.Title.Interpreter = 'latex';
C.Title.FontSize = 20;

%%% Labels and stuff.
xlim([0 1.01])
ylim([0.75 1.5])
xlabel('VW Block Length ($L_{VW} / h^*_F$)', 'Interpreter', 'latex')
ylabel('Fault Averaged $(a/b)_M$', 'Interpreter', 'latex')
ax = gca;
% ax.FontSize = 18;
ax.TickLabelInterpreter = 'latex';
ax.XTick = 0:0.1:1;
ab_M1 = ax.YLim(1);
ab_M2 = ax.YLim(2);

%%% Add eta axis on the right.
yyaxis right
b = 0.01;
a_VW = 0.0069;
a_VS = 0.0154;
ab_VW = a_VW/b;
ab_VS = a_VS/b;
eta1 = (ab_M1 - ab_VS)/(ab_VW - ab_VS);
eta2 = (ab_M2 - ab_VS)/(ab_VW - ab_VS);
ylim([eta2 eta1])

ax2 = gca;
ax2.YDir = 'reverse';
ax2.YColor = 'k';
ax2.YLabel.Interpreter = 'latex';
ax2.YLabel.String = 'VW Fraction of Fault ($\eta$)';