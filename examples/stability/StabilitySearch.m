function [R, p, L_star] = StabilitySearch(ab, beta, mu_0, geom, Burial_Star, tol, calc, R)

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
        L_star = (a + b)/2;
   
    %%% Set up the model.
        if exist('R', 'var')
            [R, p] = StabilitySetup(ab, L_star, beta, mu_0, geom, Burial_Star, R);
        else
            [R, p] = StabilitySetup(ab, L_star, beta, mu_0, geom, Burial_Star);
        end
    
    %%% Conduct the stability calculations.
        Out_S = StabilityAnalysis(p, calc, geom);
        stab = Out_S.Stability;

    %%% Update the search interval.       
        if stab == stab_a
            a = L_star;
            % stab_a = stab;
        else
            b = L_star;
        end

    %%% Increment the step counter.
        q = q + 1;
    end

end