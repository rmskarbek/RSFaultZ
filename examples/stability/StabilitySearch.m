function [R, p, L_hat, tol] = StabilitySearch(ab, beta, mu_0, geom, Burial_hat, calc, type, R)

%%% For a given set of RSF and elastic parameters, this function carries out a search
%%% for the critical fault length where the stability changes. This function 
%%% implements a bisection method to determine the critical fault length to within
%%% some arbitrary tolerance.

%%% For fault length normalized by the full-space critical wavelength L_star =
%%% L/h_F_star, assume that stability changes in the interval [a, b] = [0, 2].
    a = 0;
    b = 2;

    stab_a = 1;     % System is stable
    % stab_b = 0;     % System is unstable

%%% The smallest allowable search interval is the grid spacing. So set up a grid to get this
%%% value. Use a dummy value for L_hat.
    L_hat1 = 1;
    if exist('R', 'var')
        [R, p] = StabilitySetup(ab, L_hat1, beta, mu_0, geom, Burial_hat, type, R);
    else
        [R, p] = StabilitySetup(ab, L_hat1, beta, mu_0, geom, Burial_hat, type);
    end
%%% Search tolerance is defined as 1/2 of the search interval. Account for the different
%%% scaling for thin layer geometry.
    if strcmp(geom, 'Layer') == 1
        d = Burial_hat*p.Friction.L_b;
        H_s = (8*pi*d*p.Friction.H_star)^(1/2);
        tol = 0.5*(p.Geometry.GridSpacing/H_s);
    else
        tol = 0.5*(p.Geometry.GridSpacing/p.Friction.H_star);
    end        

%%% Maximum number of steps.
    q_max = 100;
    q = 1;
    while q < q_max

    %%% Stop the loop if the half the length of the current interval is less than the
    %%% tolerance.
        if (b - a)/2 < tol
            break
        end

    %%% Check the stability at the midpoint of current interval.
        L_hat = (a + b)/2;
   
    %%% Set up the model.
        [R, p] = StabilitySetup(ab, L_hat, beta, mu_0, geom, Burial_hat, type, R);
    
    %%% Conduct the stability calculations.
        Out_S = StabilityAnalysis(p, calc, geom);
        stab = Out_S.Stability;

    %%% Update the search interval.       
        if stab == stab_a
            a = L_hat;
            % stab_a = stab;
        else
            b = L_hat;
        end

    %%% Increment the step counter.
        q = q + 1;
    end

end