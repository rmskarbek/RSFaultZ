%%%---------------------------------------------------------------------%%%

%%% This script executes two calculations:

%%% Calculation 1: Numerically determine the critical fault length, by 
%%% finding the smallest fault length that is unstable. Do the calculation
%%% for the full-space and thin layer geometries.

%%% Calculation 2: For the thin layer geometry, determine the wavelength of
%%% the 1st unstable mode as a function of fault length. This wavelength 
%%% approaches H_s_L as the fault length increases.

%%%---------------------------------------------------------------------%%%

%%% This script executes the calculations according to the following steps:
%%% 1. Sets up ranges of parameter values in this script.

%%% 2. Calls StabilityRun1.m: executes the stability analysis for each set
%%% of parameter values; this occurs in two steps.

%%%     2a. StabilitySetup.m: Creates an instance of RSFaultZ and sets the
%%%     parameter values.

%%%     2b. StabilityAnalysis.m: Executes the numerical stability analysis
%%%     using the parameters structure that is created by StabilitySetup.m

%%% 3. Generates Figure 2 in Skarbek, 2024 JGR.

%%%---------------------------------------------------------------------%%%

%%% This script takes ~3 minutes to execute on Rob's computer: 3.7 GHz 
%%% processor, 64 GB RAM, Linux OS. Essentially all of the execution time
%%% is during Calculation 2 because it requires full eigenvalues; whereas
%%% Calculation 1 only requires finding the eigenvalue with the largest
%%% real part. 

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Calculation 1: Numerically determine the critical fault length, by 
%%% finding the smallest fault length that is unstable. Do the calculation
%%% for the full-space and thin layer geometries.

%%% Neither geometry has a dip angle. The value of mu_0 does not affect
%%% these calculations, but needs to be set to some value.
beta = nan;
mu_0 = 0.6;

%%% Values of a/b.
ab1 = (0.1:0.1:0.9)';

%%% Run the stability analysis. This calculation requires computing only
%%% the eigenvalue with the largest real part.
calc = 'boundary';

%%% Full space. Values of fault length as Length_star = FaultLength/h_star_F.
geom = 'Full Space';
Burial_Star = nan;
Length_starF = (0.35:0.002:0.368)';
[T_F, P_F] = StabilityRun1(ab1, Length_starF, calc, geom, beta, mu_0, Burial_Star);

%%% Thin layer. Values of fault length as Length_star = FaultLength/h_star_L,
%%% and burial depth as Burial_Star = BurialDepth/h_star_L.
geom = 'Layer';
Burial_Star = 0.01;
Length_starL = (0.49:0.002:0.508)';
[T_L, P_L] = StabilityRun1(ab1, Length_starL, calc, geom, beta, mu_0, Burial_Star);

%%% This loop pulls out the desired fault length for each value of a/b.
L_star_F = nan(numel(ab1),1);
L_star_L = nan(numel(ab1),1);
for i = 1:numel(ab1)
    
    for j = 1:numel(Length_starF)
        if T_F.(j){i,1}.Stability == 0
            L_star_F(i,1) = P_F.(j){i,1}.Geometry.FaultLength;
            break
        end
    end

    for j = 1:numel(Length_starL)
        if T_L.(j){i,1}.Stability == 0
            L_star_L(i,1) = P_L.(j){i,1}.Geometry.FaultLength;
            break
        end
    end
end

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Calculation 2: For the thin layer geometry, determine the wavelength of
%%% the 1st unstable mode as a function of fault length. This wavelength 
%%% approaches H_s_L as the fault length increases.
ab2 = 0.5;
Length_starL2 = (0.2:0.2:10)';
calc = 'full';
[T_L2, P_L2] = StabilityRun1(ab2, Length_starL2, calc, geom, beta, mu_0, Burial_Star);

%%% Analytic result for wavelengths as a function of fault length.
H_s_F = P_L2.(1){1,1}.Friction.H_star;
L_b = P_L2.(1){1,1}.Friction.L_b;
d = Burial_Star*L_b;
H_s_L = (8*pi*d.*H_s_F).^(1/2);

L_star = linspace(0.5, max(Length_starL2), 500)';
Modes = floor(2*L_star);
i_zero = Modes == 0;
Modes(i_zero) = nan;
WaveLength = 2*L_star*H_s_L./Modes;

%%% Wavelengths from numerical calculation.
WaveLength3 = nan(numel(Length_starL2),1);
j = 1;
for i = 1:numel(Length_starL2)    
    if T_L2.(i){j,1}.Stability == 0
        [~, kk_pos] = findpeaks(real(T_L2.(i){j,1}.vector_1));
        n_pos = numel(kk_pos);
        [~, kk_neg] = findpeaks(-real(T_L2.(i){j,1}.vector_1));
        n_neg = numel(kk_neg);
        
        if n_pos == 1 && n_neg == 0 || n_pos == 0 && n_neg == 1
            WaveLength3(i,1) = 2*P_L2.(i){j,1}.Geometry.FaultLength;
        elseif n_pos == 1 && n_neg == 1
            Xi = P_L2.(i){j,1}.Geometry.FaultCoordinates;
            WaveLength3(i,1) = 2*abs(Xi(kk_pos) - Xi(kk_neg));
        elseif n_pos == 1 && n_neg == 2
            Xi = P_L2.(i){j,1}.Geometry.FaultCoordinates;
            WaveLength3(i,1) = diff(Xi(kk_neg));
        else
            WaveLength3(i,1) = mean(T_L2.(i){j,1}.WaveLengths);
        end
    end
end

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Analytic values of h_star_L. Use any p structure for elastic
%%% parameters.
p = P_F.(1){1,1};
ShearMod = p.Material.ShearModulus;
Poisson = p.Material.PoissonRatio;
M = ShearMod/(1 - Poisson);
EffStress_0 = p.Material.BackgroundStress;

b = p.Friction.b;
a = b*linspace(0.01, 0.95, 100)';
d_c = p.Friction.d_c;

L_b = p.Friction.L_b;
H_s_F = pi*M*d_c./(EffStress_0*(b-a));

d = Burial_Star*L_b;
H_s_L1 = (8*pi*d.*H_s_F).^(1/2);
L_bh = (2*M*d*d_c/(b*EffStress_0)).^(1/2);

%%% Plot the results from Calculation 1.
figure;
tiledlayout(1,2);
nexttile;

plot(a/b, 1e-3*H_s_L1/2, 'c', 'LineWidth', 3);
hold on;
plot(a/b, 1e-3*exp(-1)*H_s_F, 'm', 'LineWidth', 3);
plot(ab1, 1e-3*L_star_L, 'ks', 'MarkerFaceColor', 'c', 'MarkerSize', 15)
plot(ab1, 1e-3*L_star_F, 'ks', 'MarkerFaceColor', 'm', 'MarkerSize', 15)
hold off

%%% Labels and such.
xlabel('$a/b$', 'Interpreter', 'latex')
ylabel('Critical Fault Length $L^*$ (km)', 'Interpreter', 'latex')

xlim([0 0.95])
ylim([0 35])

ax = gca;
ax.FontSize = 20;
ax.XTick = 0:0.1:0.9;
ax.YTick = 0:5:35;
ax.TickLabelInterpreter = 'latex';

lgd = legend('$L^*_L = h^*_L/2$', '$L^*_F = h^*_F/e$', 'Location',...
    'northwest', 'Interpreter', 'latex');
text(nexttile(1), 0.9, 0.9, 'A', 'units', 'normalized', 'FontSize', 30,...
    'Interpreter', 'latex')

%%%---------------------------------------------------------------------%%%
%%% Plot the results from Calculation 2.
nexttile;
plot(L_star, WaveLength/H_s_L, 'k', 'LineWidth', 2)
hold on
plot([0 max(L_star)], [1 1], 'k--', 'LineWidth', 2)
plot(Length_starL2, WaveLength3/H_s_L, 'ks', 'MarkerFaceColor', 'c',...
    'MarkerSize', 10)
plot([0.5 0.5], [0.5 2.5], 'k--', 'LineWidth', 2)
hold off

%%% Labels and such.
xlabel('Fault Length ($L / h_L^*$)', 'Interpreter', 'latex')
ylabel('Wavelength ($\lambda_{nT} / h_L^*$)', 'Interpreter', 'latex')

xlim([0 10])
ylim([0.5 2.5])

ax = gca;
ax.FontSize = 20;
ax.XTick = [0, 0.5, 1:max(L_star)];
ax.YTick = 0.5:0.5:2.5;
ax.TickLabelInterpreter = 'latex';

text(nexttile(2), 0.9, 0.9, 'B', 'units', 'normalized', 'FontSize', 30,...
    'Interpreter', 'latex')