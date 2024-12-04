function [T, P] = BlocksRun(aminusb_VW, aminusb_VS, ab_M, FaultLength, Lhat_VW, geom,...
    calc, type, effstressF, symmetry)

%%% This function sets up and runs different calculations for heterogeneous block faults 
%%% using RSFaultZ.

%%% INPUT -----------------------------------------------------------------------------%%%
% aminusb_VW  - A singe value for RSF parameter (a-b) for velocity-weakening material.

% aminusb_VM  - A singe value for RSF parameter (a-b) for velocity-strengthening material.

% ab_M        - A vector of values for the fault's mean value of (a/b).

% FaultLength - The target length of the fault in meters.

% Lhat_VW     - A vector of values for the length of velocity-weakening blocks, normalized
%               by h*_F.

% geom        - A string to determine the geometry of the fault system. Possible values
%               are from RSFaultZ. However, here only geom = 'Thrust Fault' is used.

% calc        - A string that is used for stability calculations. This is input into
%               StabilityAnalysis.m, see that code for more information. This variable
%               must have a value, but is only used if type = 'stability'.

% type        - A string that determines what type of calculation to run. Possible values:
%               type = 'simulation' - run a simulation of the nonlinear governing
%                       equations.
%               type = 'stability' - conduct a numerical linear stability analysis.
%               type = 'grid' - just generate the numerical grid.

% effstressF  - A string that determines how to hand the effective stress on the fault.
%               effstressF = 'constant' - effective stress everywhere on the dipping fault
%                             is set equal to the lithostatic stress at the midpoint of 
%                             the fault.
%               effstressF = 'variable' - the effective stress on the dipping fault varies
%                             with depth accordind to a pore pressure ratio of 0.75;

% symmetry    - A string that determines the symmetry of the heterogeneous block grid.
%               symmetry = 'asymmetric_VW' - VW block at up-dip edge and VS block at
%                           down-dip edge.
%               symmetry = 'asymmetric_VS' - VS block at up-dip edge and VW block at
%                           down-dip edge.
%               symmetry = 'asymmetric_VW' - VW blocks at up-dip and down-dip edges.
%%%------------------------------------------------------------------------------------%%%


%%% OUTPUT ----------------------------------------------------------------------------%%%
% T         - A table that stores the results of the desired type of calculation for each
%             pair of values given by the vectors ab_M and Lhat_VW. When type =
%             'simulation', ony reduced simulation results are stored in T. The full
%             results are automatically saved into a .mat file.

% p         - A table that stores an RSFaultZ parameters structure for each pair of values
%             given by the vectors ab_M and Lhat_VW.
%%%------------------------------------------------------------------------------------%%%

%%%---------------------------------------------------------------------%%%
%%% Create tables to store the results. Access elements like: T{j,i}{1,1}.
%%%---------------------------------------------------------------------%%%
%%% Values of (a/b)_M.
N_aM = numel(ab_M);

%%% Values of Lhat_VW = Length_VW/L_star.
N_L = numel(Lhat_VW);

%%% Create names for the different variables.
ab_MNames = cell(N_aM, 1);
% ab_WNames = cell(N_aW, 1);
LhatNames = cell(N_L, 1);

%%% Names for (a/b)_M values.
for q = 1:N_aM
    if ab_M(q) < 1
        ab_MNames{q,1} = ['abM_0', num2str(1e3*ab_M(q))];
    else
        ab_MNames{q,1} = ['abM_', num2str(1e3*ab_M(q))];
    end
end

%%% Names for Lhat_VW values.
for q = 1:N_L
    if Lhat_VW(q) < 0.1
        LhatNames{q,1} = ['Lhat_0', num2str(1e2*Lhat_VW(q))];
    elseif mod(Lhat_VW(q), 0.1) ~= 0
        LhatNames{q,1} = ['Lhat_', num2str(1e3*round(Lhat_VW(q),3))];
    else
        LhatNames{q,1} = ['Lhat_', num2str(1e1*Lhat_VW(q))];
    end
end

rowNames = ab_MNames;
columnNames = LhatNames;

N = numel(rowNames);
M = numel(columnNames);
T = table('Size', [N M], 'VariableTypes', repmat({'cell'}, 1, M),...
    'VariableNames', columnNames);

T.Properties.RowNames = rowNames;

%%% Create a table to store parameter structures.
P = T;

%%%---------------------------------------------------------------------%%%
%%% Run the calculations. Loop i over columns, and j over rows.
%%%---------------------------------------------------------------------%%%
k = 1;
y = 1;
q = 1;
for i = 1:M
    for j = 1:N

    %%% Determine which variables to loop over.
        switch rowNames{1,1}(1:3)
            case 'abM' % ab_MNames
                k = j;
            case 'Lha' % LhatNames
                q = j;            
        end
    
        switch columnNames{1,1}(1:3)
            case 'abM' % ab_MNames
                k = i;
            case 'Lha' % LhatNames
                q = i;
        end

        switch type
            case 'grid'

                G = BlocksSetup(aminusb_VW(y), aminusb_VS, ab_M(k), Lhat_VW(q),...
                    FaultLength, geom, effstressF, symmetry, type);
                T{j,i} = {G.FaultLength};
                P{j,i} = {G};

            case 'simulation'

                if exist('R', 'var')
                    [R, p] = BlocksSetup(aminusb_VW(y), aminusb_VS, ab_M(k), Lhat_VW(q),...
                        FaultLength, geom, effstressF, symmetry, type, R);
                else
                    [R, p] = BlocksSetup(aminusb_VW(y), aminusb_VS, ab_M(k), Lhat_VW(q),...
                        FaultLength, geom, effstressF, symmetry, type);
                end

            %%% Run the simulation.
                R.StartButton.ButtonPushedFcn([],[]);
            
            %%% Downsample factor.
                ds = 20;
                switch p.Options.DataManagement
                    case 'ode'
                    %%% Access the simulation output.
                        SimData = struct('p', p, 'Output', R.Out);
                    %%% Maximum slip velocity in VW regions.
                        MaxVW = max(max(SimData.Output.Vel(:,p.Geometry.i_VW), [], 2));
                    %%% Maximum slip velocity in VW centers.
                        MaxVWc = max(max(SimData.Output.Vel(:,p.Geometry.i_VW_center), [], 2));
                    %%% Maximum slip velocity in center of all VS regions.
                        MaxVSc = max(max(SimData.Output.Vel(:,p.Geometry.i_VS_center), [], 2));
                    %%% Downsample the output.                        
                        SimData.Output.Time = downsample(SimData.Output.Time, ds);
                        SimData.Output.Vel = downsample(SimData.Output.Vel, ds);
                        SimData.Output.ShearStress = downsample(SimData.Output.ShearStress, ds);
                        SimData.Output.NormalStress = downsample(SimData.Output.NormalStress, ds);
                        SimData.Output.State = downsample(SimData.Output.State, ds);
                        SimData.Output.Slip = downsample(SimData.Output.Slip, ds);
                    %%% Save the simulation data and then clear it to free up memory.                
                        save([ab_MNames{k}, '_', LhatNames{q}], 'SimData', 'MaxVW', 'MaxVWc',...
                            'MaxVSc');
                        clear SimData                
                
                    case 'dat'
                        Vel = SaveStreamData('Read', 'Vel.dat')';
                    %%% Maximum slip velocity in VW regions.
                        MaxVW = max(max(Vel(:,p.Geometry.i_VW), [], 2));
                    %%% Maximum slip velocity in VW centers.
                        MaxVWc = max(max(Vel(:,p.Geometry.i_VW_center), [], 2));
                    %%% Maximum slip velocity in center of all VS regions.
                        MaxVSc = max(max(Vel(:,p.Geometry.i_VS_center), [], 2));
        
                        Vel = downsample(Vel, ds);
                        Time = downsample(SaveStreamData('Read', 'Time.dat')', ds);
                        NormalStress = downsample(SaveStreamData('Read', 'Sigma.dat')', ds);
                        ShearStress = downsample(SaveStreamData('Read', 'Tau.dat')', ds);
                        State = downsample(SaveStreamData('Read', 'State.dat')', ds);
        
                        SimData.p = p;
                        Output = struct('NormalStress', NormalStress, 'ShearStress', ShearStress,...
                            'State', State, 'Vel', Vel, 'Time', Time);
                        clear Vel Time NormalStress ShearStress State
                        SimData.Output = Output;
                        clear Output
                        save([ab_MNames{k}, '_', LhatNames{q}], 'SimData', 'MaxVW', 'MaxVWc',...
                            'MaxVSc');
                        clear SimData
                    %%% Clean up the .dat files.
                        delete Sigma.dat Tau.dat Time.dat Vel.dat State.dat
                end

            %%% Assemble.
                V = struct('MaxVW', MaxVW, 'MaxVWc', MaxVWc, 'MaxVSc', MaxVSc);
                T{j,i} = {V};
                P{j,i} = {p};

            case 'stability'

             %%% Conduct the stability search for the critical VW block length.                
            %     if exist('R', 'var')
            %         [R, p, L_star_VW, tol] = P_StabilitySearch(aminusb_VW(y), aminusb_VS,...
            %             ab_M(k), FaultLength, geom, effstressF, symmetry, calc, R);
            %     else
            %         [R, p, L_star_VW, tol] = P_StabilitySearch(aminusb_VW(y), aminusb_VS,...
            %             ab_M(k), FaultLength, geom, effstressF, symmetry, calc);
            %     end
            % %%% Tack the tolerance onto p.
            %     p.tol = tol;
            % 
            % %%% Store the critical VW block length directly from the RSFaultZ parameters
            % %%% structure, and store the structure as well.
            %     T{j,i} = {L_star_VW};
            %     P{j,i} = {p};

            %%% Locate the p-stability boundary. Assume stability = 1 with Lhat_VW = 0.2.
                tol = 0.005;
                qq = 0;
                Lhat_VW = 0.2;
                dLhat_VW = 0.2;
                stab_0 = 1;
                while abs(dLhat_VW) > tol
                    qq = qq + 1;

                %%% Set up the model.
                    if exist('R', 'var')
                        [R, p] = BlocksSetup(aminusb_VW(y), aminusb_VS, ab_M(k), Lhat_VW(q),...
                            FaultLength, geom, effstressF, symmetry, type, R);
                    else
                        [R, p] = BlocksSetup(aminusb_VW(y), aminusb_VS, ab_M(k), Lhat_VW(q),...
                            FaultLength, geom, effstressF, symmetry, type);
                    end

                %%% Conduct the stability calculations.
                    Out_S = StabilityAnalysis(p, calc, geom);
                    stab = Out_S.Stability;

                %%% Check is the stability changed and update the fault length
                %%% accordingly.                    
                    stab_change = stab ~= stab_0;
                    stab_0 = stab;

                    if stab_change == 1
                    %%% Stability changed, so reduce the step size and change directions.
                        if stab == 1
                            dLhat_VW = abs(dLhat_VW/2);
                            Lhat_VW = Lhat_VW + dLhat_VW;

                        else
                            dLhat_VW = -abs(dLhat_VW/2);
                            Lhat_VW = Lhat_VW + dLhat_VW;
                        end
                    else
                    %%% Stability did not change, so keep moving in the same direction
                    %%% with the same step size.
                        if stab == 1
                            dLhat_VW = abs(dLhat_VW);
                            Lhat_VW = Lhat_VW + dLhat_VW;    
                        else
                            dLhat_VW = -abs(dLhat_VW);
                            Lhat_VW = Lhat_VW + dLhat_VW;
                        end
                    end
                end

            %%% This puts Lhat_VW known to within p/m tol.
                T{j,i} = {Lhat_VW};
                P{j,i} = {p};
        end
    end
end

%%% Close the RSFaultZ instance.
if exist('R', 'var')
    delete(R);
end