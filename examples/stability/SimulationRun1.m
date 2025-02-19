function [T, P] = SimulationRun1(ab, Length_hat, geom, beta, mu_0, calc,...
    Burial_hat, EffStress_0)
%%% This function is altered from StabilityRun1.m, to run full simulations.

%%% Values of a/b.
N_a = numel(ab);

%%% Values of Length_star = FaultLength/h_star.
N_L = numel(Length_hat);

%%% Create names for the different variables.
abNames = cell(N_a, 1);
LhatNames = cell(N_L, 1);

%%% Names for a/b values.
for q = 1:N_a
    abNames{q,1} = ['ab_0', num2str(10*ab(q))];
end

%%% Names for Length_star values.
for q = 1:N_L
    LhatNames{q,1} = ['Lhat_', num2str(Length_hat(q))];
end

%%%---------------------------------------------------------------------%%%
%%%---------------------------------------------------------------------%%%

%%% Each calculation varies two parameters out of: a/b, L_star, beta, mu_0.
%%% For full space and layer geometries, only a/b and L_star are varied,
%%% since the values of beta and mu_0 do not effect the calculation.

%%% Create a table to store the results. Access elements like: T{j,i}{1,1}.
switch geom

    case {'Full Space', 'Layer'}
        rowNames = abNames;
        columnNames = LhatNames;
        
        N = numel(rowNames);
        M = numel(columnNames);
        T = table('Size', [N M], 'VariableTypes', repmat({'cell'}, 1, M),...
            'VariableNames', columnNames);

    case {'Thrust Fault', 'Normal Fault'}
    %%% This is for the L_min calculation for half space.
        rowNames = betaNames;
        columnNames = BstarNames;

        N = numel(rowNames);
        % M = 1;
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
            case 'L' % LHatNames
                q = j;
            case 'b' % betaNames
                w = j;
            case 'm' % MuNames
                y = j;
            case 'B' % BhatNames
                z = j;
        end
    
        switch columnNames{1,1}(1)
            case 'a' % abNames
                k = i;
            case 'L' % LHatNames
                q = i;
            case 'b' % betaNames
                w = i;
            case 'm' % MuNames
                y = i;
            case 'B' % BhatNames
                z = i;
        end

        switch calc
            case 'length'

            %%% Setup RSFaultZ.
                if exist('R', 'var')
                    [R, p] = SimulationSetup(ab(k), Length_hat(q), beta(w),...
                        mu_0(y), geom, Burial_hat(z), EffStress_0, R);
                else
                    [R, p] = SimulationSetup(ab(k), Length_hat(q), beta(w),...
                        mu_0(y), geom, Burial_hat(z), EffStress_0);
                end

            %%% Alter the initial velocityconditions to control the perturbation
            %%% length. Here, only a singe node at the center of the fault is
            %%% perturbed.
                N_g = p.Geometry.GridPoints;
                N_p = round(N_g/2);
                v_plate = p.Material.PlateRate;

                vel_i = R.p.Geometry.InitConditions(1:N_g,1);
            %%% Perturbation at the center of the fault.
                % vel_i(1:N_p-1) = v_plate;
                % vel_i(N_p+1:end) = v_plate;
                % vel_i(N_p) = 0.99*v_plate;
            %%% Perturbation at the edge of the fault.
                vel_i(1:end) = v_plate;
                vel_i(end) = 0.99*v_plate;
                R.p.Geometry.InitConditions(1:N_g,1) = vel_i;

            %%% Run the simulation.
                R.StartButton.ButtonPushedFcn([],[]);

            %%% Access the simulation output.
                SimData = struct('p', R.p, 'Output', R.Out);
                Time = downsample(SimData.Output.Time, 10);
                V_max = downsample(max(SimData.Output.Vel,[],2), 10);             
                V = struct('Time', Time, 'Vmax', V_max);
                T{j,i} = {V};
                P{j,i} = {R.p};

            case 'perturbation'

%%% Setup RSFaultZ. Use a constant fault length, Length_star = 1. Input
%%% values of Length_star will be used to control the perturbation length.
                Length_starP = 1;
                if exist('R', 'var')
                    [R, p] = StabilitySetup(ab(k), Length_starP, beta(w),...
                        mu_0(y), geom, Burial_hat(z), R);
                else
                    [R, p] = StabilitySetup(ab(k), Length_starP, beta(w),...
                        mu_0(y), geom, Burial_hat(z));
                end
 
%%% Alter the initial state variable conditions to control the perturbation
%%% length.
                N_g = p.Geometry.GridPoints;
                N_p = round(Length_hat(q)*N_g);
                N_steady = round((N_g - N_p)/2);                
                v_plate = p.Material.PlateRate;

                vel_i = R.p.Geometry.InitConditions(1:N_g,1);
                vel_i(1:N_steady) = v_plate;
                vel_i(N_steady+N_p+1:end) = v_plate;
                R.p.Geometry.InitConditions(1:N_g,1) = vel_i;

%%% Run the simulation.
                R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
                SimData = struct('p', R.p, 'Output', R.Out);
                Time = downsample(SimData.Output.Time, 10);
                V_max = downsample(max(SimData.Output.Vel,[],2), 10);             
                V = struct('Time', Time, 'Vmax', V_max);
                T{j,i} = {V};
                P{j,i} = {R.p};
        end
    end
end

%%% Close the RSFaultZ instance.
delete(R);