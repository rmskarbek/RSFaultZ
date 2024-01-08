function [Beta, Shear, G_thrust, G_normal] = StressingRates

%%% This script computes normalized stressing rates on a dipping fault with
%%% a uniform velocity distribution.

%%% This script executes in 0.002 s on Rob's computer: 3.7 GHz processor,
%%% 64 GB RAM, Linux OS; and carries out calculations for Figure 5B in 
%%% Skarbek, 2024 JGR.

%%%---------------------------------------------------------------------%%%

%%% Dip angles [deg].
    Beta = linspace(0, 90, 100)';

%%% Fault Length [m].
    L = 100;

%%% Surface slope [rad].
    alpha = 0;

%%% Burial Depth [m].
    d = 0;

%%% Define grid, although only one fault element is required for this
%%% calculation, two are used for proper input into WedgeStress_HODLR.m
    N = 2;
    dxi = L;
    Xi = 0*dxi + (dxi/2: dxi: N*dxi - dxi/2)';

%%% Eavluate the stress change functions for just the first element.
    i = 1;
    j = 1;
    Shear = nan(numel(Beta),1);
    Normal = nan(numel(Beta),1);

    for k = 1:numel(Beta)
        S = WedgeStress_HODLR(i, j, Xi', Beta(k), alpha, d, 1);
        N = WedgeStress_HODLR(i, j, Xi', Beta(k), alpha, d, 0);

        Shear(k,1) = L*S;
        Normal(k,1) = L*N;
    end
      
%%% Cumpute the stressing rates.
    mu_0 = [0.2, 0.6, 1];
    G_thrust = abs(Shear - mu_0.*Normal);
    G_normal = abs(Shear + mu_0.*Normal);