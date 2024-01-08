function [Mu, DH, DH_max, AB, AB_max, K_c, K_max] = F6_LayerVS

%%%---------------------------------------------------------------------%%%
%%% This function carries out linear stability calculations for an
%%% infinitely long, velocity-strengthening fault that is parallel to a
%%% traction-free surface. 

%%% This function will automatically generate Figure 6 in Skarbek, 2024,
%%% JGR.

%%% This function takes ~50 s to run on Rob's computer: 3.7 GHz processor,
%%% 64 GB RAM, Linux OS.
%%%---------------------------------------------------------------------%%%

%%% Normalized burial depth: d_hat = d/L_b;
DH = logspace(0, 4, 200)';

%%% Values for mu_0 and dimensionless wavenumber k_hat = k*d for the 
%%% calculation. Numerical tests indicate that when there is an unstable VS
%%% region, it occurs for values of k of around 1 (Figure B1).
Mu = linspace(0,1,200)';
k_hat = linspace(0,10,100)';

%%% Run the calculation over the values of DH, and Mu. For each value of
%%% DH and Mu, the loop will find the maximum value of a/b > 1 that yields
%%% an eigenvalue with a postive real part. That value of a/b is then the
%%% critical VS value for that (d_hat, mu_0) pair.
AB = nan(numel(Mu), numel(DH));

%%% Also store the wavenumber at the critical value of a/b. For 1 < a/b <
%%% (a/b)_c, there will be a range of unstable wavenumbers. For a/b =
%%% (a/b)_c there is only one unstable wavenumber.
K_c = nan(numel(Mu), numel(DH));

%%% Execute the calculation loop.
for j = 1:numel(DH)
    d_hat = DH(j);

%%% Stress change functions in Fourier space. Only depend on d_hat and
%%% k_hat.
    T_v = -(k_hat/(2*d_hat)).*(1 - exp(-2*k_hat).*(1 - 2*k_hat + 2*k_hat.^2));
    N_v = -1i*(k_hat.^3/d_hat).*exp(-2*k_hat);

    for i = 1:numel(Mu)
        mu_0 = Mu(i);
%%% Initial value of a/b is VW, so that it will produce an unstable
%%% eigenvalue.
        ab = 0.99999;
    
%%% Compute the eigenvalues, the positive solution always has a larger real
%%% part, so just use that one. The negative solution is commented out as
%%% lambda_2.
        Gamma = T_v - mu_0*N_v;
        B = 1 - (1./ab).*(Gamma + 1);
        C = -(1./ab).*Gamma;
           
        lambda_1 = (1/2)*(-B + sqrt(B.^2 - 4*C));
        % lambda_2 = (1/2)*(-B - sqrt(B.^2 - 4*C));
        
%%% Now increase values of a/b until the real part of lambda_1 becomes
%%% negative.
        q = 0;
        while max(real(lambda_1)) > 0
            q = q + 1;
            ab = ab + 0.00001;
        
            B = 1 - (1./ab).*(Gamma + 1);
            C = -(1./ab).*Gamma;
            lambda_1 = (1/2)*(-B + sqrt(B.^2 - 4*C));
        end

%%% The last value of q was stable, so record the second to last value, as
%%% the largest value of a/b > 1 that is unstable.
        AB(i,j) = ab - 0.00001;

%%% Compute lambda_1 again to find the critical wavelength.
        B = 1 - (1./ab).*(Gamma + 1);
        C = -(1./ab).*Gamma;
        lambda_1 = (1/2)*(-B + sqrt(B.^2 - 4*C));
        [~, ii] = findpeaks(real(lambda_1));
        if isempty(ii) == 0
            K_c(i,j) = k_hat(ii);
        end
    end
end

%%% For each value of mu_0, find the value of d_hat that corresponds to the
%%% maximum value of (a/b)_c.
[AB_max, I] = max(AB,[],2);
DH_max = DH(I);
K_max = K_c(sub2ind(size(K_c), (1:200)', I));

%%% Analytic equations for minimum depth, and depth of max (a/b)_c.
% Mu1 = -((2*exp(1) - 1)./(1 + 2*(4*DH-1)*exp(1))).^(1/2)*(1-2*exp(1));
% Mu2 = -((exp(2) - 1)./(1 + (2*DH-1)*exp(2))).^(1/2)*(1-exp(2));

%%% Approximate analytic expression for the minimum depth.
D_min = (2./Mu).^2;

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Contour plot of (a/b)_c results. Remove velocity-weakening values.
ii = AB == 0.99999;
AB(ii) = nan;
figure;
tiledlayout(1,2);
nexttile
% levels = linspace(1, max(max(AB)), 12)';
[~, C2] = contourf(repmat(DH',numel(Mu),1), repmat(Mu,1,numel(DH)), AB);
hold on
plot(D_min, Mu, 'k--', 'LineWidth', 3);
hold off
xlim([1 1e4])

C = colorbar('TickLabelInterpreter', 'latex');
colormap(cool)
C.Ticks = C2.LevelList;
C.TickLabels(2:2:end) = {''};
% C.Ticks = [1.00001, C.Ticks];
% C.TickLabels{1} = '1';

C.Title.String = '$(a/b)_c$';
C.Title.Interpreter = 'latex';
C.Title.FontSize = 24;

% hold on
% plot(LH, Mu_min, 'c')
% plot(LH_max, Mu, 'r')
% hold off

lgd = legend('', '$\mu_0 = (4 L_b / d)^{1/2}$', 'Location', 'southwest',...
    'Interpreter', 'latex');

ax = gca;
ax.XScale = 'log';
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';

xlabel('Normalized Burial Depth ($d / L_b$)', 'Interpreter', 'latex')
ylabel('Friction Coefficient ($\mu_0$)', 'Interpreter', 'latex')

text(0.05, 0.93, 'A', 'units', 'normalized', 'FontSize', 30, 'Interpreter',...
    'latex')

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%
%%% Contour plot of critical wavelength values.
nexttile
contourf(repmat(DH',numel(Mu),1), repmat(Mu,1,numel(DH)), K_c)
xlim([1 1e4])

C = colorbar('TickLabelInterpreter', 'latex');
colormap(cool)
C.Ticks = [0.506, C.Ticks];

C.Title.String = '$\hat{k} = k_c/d$';
C.Title.Interpreter = 'latex';
C.Title.FontSize = 24;

ax = gca;
ax.XScale = 'log';
ax.FontSize = 20;
ax.TickLabelInterpreter = 'latex';

xlabel('Normalized Burial Depth ($d / L_b$)', 'Interpreter', 'latex')
ylabel('Friction Coefficient ($\mu_0$)', 'Interpreter', 'latex')

text(0.05, 0.93, 'B', 'units', 'normalized', 'FontSize', 30, 'Interpreter',...
    'latex')