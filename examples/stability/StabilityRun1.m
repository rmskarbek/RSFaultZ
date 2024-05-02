function [T, P] = StabilityRun1(ab, Length_star, calc, geom, beta, mu_0,...
    Burial_Star, type)
%%% This function sets up and runs numerical stability analysis for
%%% different sets of input parameters. It is called by...

%%% Each call to this functoion varies two parameters out of: a/b, L_star, 
%%% beta, mu_0.

%%% Values of a/b.
N_a = numel(ab);

%%% Values of Length_star = FaultLength/h_star.
N_L = numel(Length_star);

%%% Values of dip angle.
N_b = numel(beta);

%%% Values of friction coefficient.
N_m = numel(mu_0);

%%% Values of Burial_Star.
N_B = numel(Burial_Star);

%%% Create names for the different variables.
abNames = cell(N_a, 1);
LstarNames = cell(N_L, 1);
betaNames = cell(N_b, 1);
muNames = cell(N_m, 1);
BstarNames = cell(N_B, 1);

%%% Names for a/b values.
for q = 1:N_a
    abNames{q,1} = ['ab_0', num2str(10*ab(q))];
end

%%% Names for Length_star values.
for q = 1:N_L
    LstarNames{q,1} = ['Lstar_', num2str(Length_star(q))];
end

%%% Names for beta values.
for q = 1:N_b
    betaNames{q,1} = ['beta_', num2str(beta(q))];
end

%%% Names for mu_0 values.
for q = 1:N_m
    if mu_0(q) < 1
        muNames{q,1} = ['mu_0', num2str(10*mu_0(q))];
    else
        muNames{q,1} = ['mu_', num2str(mu_0(q))];
    end
end

%%% Names for Burial_star values.
for q = 1:N_B
    BstarNames{q,1} = ['Bstar_', num2str(Burial_Star(q))];
end

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Create a table to store the results. Access elements like: T{j,i}{1,1}.
switch geom

    case {'Full Space', 'Layer'}
        rowNames = abNames;
        switch(type)
            case 'search'
                columnNames = {'L_star'};
            case 'wavelength'
                columnNames = LstarNames;
        end
        
        N = numel(rowNames);
        M = numel(columnNames);
        T = table('Size', [N M], 'VariableTypes', repmat({'cell'}, 1, M),...
            'VariableNames', columnNames);

    case {'Thrust Fault', 'Normal Fault'}
        rowNames = betaNames;
        columnNames = BstarNames;

        N = numel(rowNames);
        M = numel(columnNames);
        T = table('Size', [N M], 'VariableTypes', repmat({'cell'}, 1, M),...
            'VariableNames', columnNames);
end

T.Properties.RowNames = rowNames;

%%% Create a table to store parameter structures.
P = T;

%%% Run the calculations. Loop i over columns, and j over rows.
k = 1;
q = 1;
w = 1;
y = 1;
z = 1;
for i = 1:M
    for j = 1:N

%%% Determine which variables to loop over.
        switch rowNames{1,1}(1)
            case 'a' % abNames
                k = j;
            case 'L' % LStarNames
                q = j;
            case 'b' % betaNames
                w = j;
            case 'm' % MuNames
                y = j;
            case 'B' % BStarNames
                z = j;
        end
    
        switch columnNames{1,1}(1)
            case 'a' % abNames
                k = i;
            case 'L' % LStarNames
                q = i;
            case 'b' % betaNames
                w = i;
            case 'm' % MuNames
                y = i;
            case 'B' % BStarNames
                z = i;
        end

        switch(type)
            case 'search'
            %%% Skip the calculation if the dip angle and the burial depth are both zero.
                if beta(w) == 0 && Burial_Star(z) == 0
                    continue
                end
        
            %%% Conduct the stability search for the critical fault length.
                tol = 0.001;
                if exist('R', 'var')
                    [R, p, L_star] = StabilitySearch(ab(k), beta(w), mu_0(y), geom,...
                        Burial_Star(z), tol, calc, R);
                else
                    [R, p, L_star] = StabilitySearch(ab(k), beta(w), mu_0(y), geom,...
                        Burial_Star(z), tol, calc);
                end
            
            %%% Store the critical fault length and the RSFaultZ parameters
            %%% structure.
                T{j,i} = {L_star};
                P{j,i} = {p};

            case 'wavelength'
            %%% Set up the fault system.
                if exist('R', 'var')
                    [R, p] = StabilitySetup2(ab(k), Length_star(q), beta(w), mu_0(y), geom,...
                        Burial_Star(z), R);
                else
                    [R, p] = StabilitySetup2(ab(k), Length_star(q), beta(w), mu_0(y), geom,...
                        Burial_Star(z));
                end
 
            %%% Conduct the stability calculations.
                Out_S = StabilityAnalysis(p, calc, geom);
                T{j,i} = {Out_S};
                P{j,i} = {p};
        end

    end
end

%%% Close the RSFaultZ instance.
delete(R);
