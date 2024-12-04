function [R, p] = BlocksSetup(aminusb_VW, aminusb_VS, ab_M, Lhat_VW, FaultLength, geom,...
    effstressF, symmetry, type, varargin)

%%% This sets up a heterogeneous grid with the goal of using constant VS and VW block 
%%% lengths, and getting to a target fault length as closely as possible. In effect
%%% controlling L_VS, L_VW, and eta at the expense of fault length. Contains options to
%%% choose a symmtric or asymmetric grid.

%----------------------------------------------------------------------------------------%
%----------------------------------------------------------------------------------------%

%%% Block symmetry.
%%% 'asymmetric_VW' - VW block first
%%% 'asymmetric_VS' - VS block first
%%% 'symmetric' - begins and ends with VS block

%%% Block friction properties. b values are from the SCEC benchmarks.
b_VW = 0.01;
% a_VW = b_VW*ab_VW;
a_VW = aminusb_VW + b_VW;
ab_VW = a_VW/b_VW;

b_VS = 0.01;
% a_VS = b_VS*ab_VS;
a_VS = aminusb_VS + b_VS;
ab_VS = a_VS/b_VS;

%%% Velocity weakening fraciton. 
eta = (ab_M - ab_VS)/(ab_VW - ab_VS);

%%% Friction coefficient.
mu_ref = 0.3;

%%% d_c tuned value[m].
% d_c = 0.004;
% d_c = 0.016;
%%% d_c upper limit from den Hartog[m].
% d_c = 0.001;
%%% d_c from SCEC benchmark [m].
d_c = 0.008;

%%% RSF reference velocity [m/s].
v_ref = 1e-6;

%%% Elastic parameters.
v_shear = 3500;
Poisson = 0.35;
%%% Shear modulus [GPa]
ShearMod = 15;
M = 1e3*ShearMod/(1 - Poisson);

%%% Constants for effective stress calculation.
g = 9.81;
beta = 8;
BurialDepth = 1e3;      % [m]
%%% Constant wedge porosity.
phi_w = 0.3;
%%% Fluid and rock densities.
rho_f = 1015;
rho_r = 2650;
%%% Normalized pore pressure ratio?
lambda = 0.75;

%%% Background effective stress scenarios. Here for the different cases, different
%%% values of effective stress at specific locations of the fault are used to set the
%%% grid spacing. These locations are determined using the target fault length, which
%%% ensures that the grid spacings will be the same for all the grids.
switch effstressF
    case 'constant'
    %%% Mean value for constant stress models. Depth of middle of fault.
        z_mid = BurialDepth + (FaultLength/2)*sind(beta);       % [m]

    %%% Vertical stress at downdip edge.
        sigma_v = g*rho_r*(1 - phi_w)*z_mid;

    %%% Hydrostatic pressure at downdip edge.
        pressure_h = g*rho_f*z_mid;

    %%% Effective stress at fault middle [MPa].
        EffStress_0 = 1e-6*(1 - lambda)*(sigma_v - pressure_h);
        % EffStress_0 = 5;

    case 'variable'
    %%% Maximum effective stress (i.e. at down-dip edge of fault).
    %%% Depth of down-dip edge.
        z_dde = BurialDepth + FaultLength*sind(beta);       % [m]
        
    %%% Vertical stress at downdip edge.
        sigma_v = g*rho_r*(1 - phi_w)*z_dde;
        
    %%% Hydrostatic pressure at downdip edge.
        pressure_h = g*rho_f*z_dde;
        
    %%% Effective stress at downdip edge [MPa].
        EffStress_0 = 1e-6*(1 - lambda)*(sigma_v - pressure_h);
end

%-----------------------------------------------------------------------------------%

%%% Friction length scales [m].
H_s_F = pi*M*d_c/(EffStress_0*(b_VW - a_VW));
L_b = d_c*M/(b_VW*EffStress_0);

%-----------------------------------------------------------------------------------%
%%% Construct the grid.
%-----------------------------------------------------------------------------------%

switch type
    case {'simulation', 'grid'}
        dxi = L_b/20;

    case 'stability'
        dxi = L_b/10;
end

%%% Velocity weakening block length [m].
L_VW = Lhat_VW*H_s_F;

switch symmetry
    case {'asymmetric_VW', 'asymmetric_VS'}
        %%% Number of elements per VW block.
        N_VW = round(L_VW/dxi);
        
        %%% Number of elements per VS block.
        N_VS = round(N_VW*(1 - eta)/eta);
        
        %%% Assemble a single VS+VW block. This defines the block as starting with 
        %%% VS material.
        A_Block = [a_VS*ones(N_VS,1); a_VW*ones(N_VW,1)];
        B_Block = [b_VS*ones(N_VS,1); b_VW*ones(N_VW,1)];
        
        %%% Find the number of blocks that gets closest to the target fault length.
        N_Block = round(FaultLength/(dxi*numel(A_Block)));
        
        %%% Make sure there is at least on VS+VW block.
        if N_Block == 0
            N_Block = 1;
        end
        
        %%% Assemble the grid by replicating the single VS+VW block.
        A = repmat(A_Block, N_Block, 1);
        B = repmat(B_Block, N_Block, 1);

%-----------------------------------------------------------------------------------%

        %%% Flip the A and B vectors to start with a VW section
        switch symmetry
            case 'asymmetric_VW'
                A = flip(A);
                B = flip(B);
        end

    case 'symmetric'
        %%% Number of VW blocks.
        N_Block = round(eta*FaultLength/L_VW);
        
        %%% Number of elements per VW block.
        N_VW = round(L_VW/dxi);
        
        %%% Number of elements per VS block.
        N_VS = round(((1 - eta)*N_Block/(eta*(N_Block + 1)))*N_VW);
        
        %%% Add a check here for reducing or increasing B_VW.
        
        %%% Assemble a single VS+VW block. This defines the block as starting with 
        %%% VS material.
        A_Block = [a_VS*ones(N_VS,1); a_VW*ones(N_VW,1)];
        B_Block = [b_VS*ones(N_VS,1); b_VW*ones(N_VW,1)];
        
        %%% Make sure there is at least on VS+VW block.
        % if N_Block == 0
        %     N_Block = 1;
        % end
        
        %%% Assemble the grid by replicating the single VS+VW block a number of times B_VW,
        %%% then add the final VS block.
        A = repmat(A_Block, N_Block, 1);
        A = [A; a_VS*ones(N_VS,1)];
        
        B = repmat(B_Block, N_Block, 1);
        B = [B; b_VS*ones(N_VS,1)];
end

%%% Final fault length and values of eta and (a/b)_M. 
N = numel(A);
FaultLength = dxi*N;
ab_M = mean(A./B);

%%% Indices of VW grid points, and eta.
i_VW = (A-B) < 0;
eta = numel(A(i_VW))/N;

%%% Final block lengths.
L_VW = dxi*N_VW;
L_VS = dxi*N_VS;
Lhat_VW = dxi*N_VW/H_s_F;
        
%%% Indices of VS grid points, and centers of VS blocks.
i_VS = i_VW == 0;
switch symmetry
    case 'asymmetric_VW'
    %%% VW first.
        i_B = (0:N_Block-1)';
        i_VS_center = N_VW + round(N_VS/2) + (N_VS + N_VW)*i_B;
        i_VW_center = round(N_VW/2) + (N_VS + N_VW)*i_B;

    case 'asymmetric_VS'
    %%% VS first.
        i_B = (0:N_Block-1)';
        i_VS_center = round(N_VS/2) + (N_VS + N_VW)*i_B;
        i_VW_center = N_VS + round(N_VW/2) + (N_VS + N_VW)*i_B;

    case 'symmetric'
        i_B_VS = (0:N_Block)';
        i_VS_center = round(N_VS/2) + (N_VS + N_VW)*i_B_VS;
        i_VW_center = N_VS + round(N_VW/2) + (N_VS + N_VW)*i_B_VS(1:end-1);
end

%%% Heterogeneous length scales.
H_s_M = pi*M*d_c/(EffStress_0*(eta*(b_VW - a_VW) + (1 - eta)*(b_VS - a_VS)));
% H_s_P = pi*M*d_c*(a_VS - eta*a_VW)/...
%   (EffStress_0*(a_VW*(b_VS - a_VS) + a_VS*(b_VW - a_VW)));

%-----------------------------------------------------------------------------------%

switch effstressF
    case 'variable'
    %%% Compute variable normal stress using the grid that RSFaultZ will generate.
    %%% Distance along-dip [m].
        Xi = 0*dxi + (dxi/2: dxi: N*dxi - dxi/2)';
    %%% Zed is the fault's depth below the upper surface [m].
        Zed = BurialDepth + Xi*(sin(pi*beta/180));
    
    %%% Vertical stress [Pa]
        sigma_v = g*rho_r*(1 - phi_w)*Zed;
    
    %%% Hydrostatic pressure [Pa].
        pressure_h = g*rho_f*Zed;
    
    %%% Effective stress [MPa].
        EffStress = 1e-6*(1 - lambda)*(sigma_v - pressure_h);
end

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
switch type
    case 'grid'
    %%% Put values into a structure for grid calculations.
        R.FaultLength = FaultLength;
        R.Lhat_VW = Lhat_VW;
        R.ab_M = ab_M;
        R.eta = eta;

    %%% Need dummy value for p.
        p = [];

    case {'simulation', 'stability'}
%%% Create an instance of the app if one is not provided.
if isempty(varargin)
    R = RSFaultZ;
else
    R = varargin{1,1};
end

%%% Set options.
%%% For full space, need to activate the function that disables the thrust
%%% geometry fields.
R.StateLawDropDown.Value = 'Aging';
% R.StateLawDropDown.Value = 'Slip';
R.PlotsOnCheckBox.Value = false;

%%% Set geometry.
R.GeometryDropDown.Value = geom;
R.FaultLengthkmEditField.Value = FaultLength/1e3;

%%% Reset the burial depth is desired.
if strcmp(geom, 'Thrust Fault') == 1 || strcmp(geom, 'Normal Fault')
    R.DipAngleEditField.Value = beta;
    % BurialDepth = 0;
    R.BurialDepthkmEditField.Value = BurialDepth/1e3;
end

%%% Set normal stress to the value used for grid spacing.
R.NormalStressMPaEditField.Value = EffStress_0;

%%% Set elastic parameters.
R.vsmsEditField.Value = v_shear;
R.PoissonRatioEditField.Value = Poisson;
R.ShearModulusGPaEditField.Value = ShearMod;

%%% Set frictional parameters
R.aEditField.Value = a_VW;
R.bEditField.Value = b_VW;
R.dcmEditField.Value = d_c;
R.v0EditField.Value = v_ref;
R.mu0EditField.Value = mu_ref;

%%% Set grid spacing.
R.dximEditField.Value = dxi;

%%% Update the Grid Controls panel.
R.dximEditField.ValueChangedFcn([],[]);

%%% Set run time.
R.RunTimeyrEditField.Value = 300;

%%% Set the properties to compute grid coordinates.
R.SetPropertiesButton.ButtonPushedFcn([],[]);

%%% At this point the grid has been set using constant properties, b_VW and
%%% EffStress_0. Now overwrite the properties that are spatially variable.

switch effstressF
    case 'variable'
    %%% Set variable effective stress.
        R.p.Material.BackgroundStress = EffStress;
    %%% Set initial conditions.
        R.p.Geometry.InitConditions(3*N+1:4*N,1) = EffStress;
    %%% Set variable L_b.
        R.p.Friction.L_b = d_c*M./(b_VW*EffStress);
end

%%% Set heterogeneous friction properties.
R.p.Friction.a = A;
R.p.Friction.b = B;

%%% Set import flag to "true", to prevent altered properties from being overwritten.
R.p.Options.Import = 1;

%%% Set data management flag.
R.p.Options.DataManagement = 'dat';

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Get parameters structure.
p = R.p;

%%% Add extra geometrical parameters.
p.Geometry.Length_VW = L_VW;
p.Geometry.Length_VS = L_VS;
p.Geometry.N_blocks = N_Block;
p.Geometry.i_VW = i_VW;
p.Geometry.i_VS = i_VS;
p.Geometry.i_VW_center = i_VW_center;
p.Geometry.i_VS_center = i_VS_center;

%%% Add extra friction parameters.
p.Friction.ab_M = ab_M;
p.Friction.eta = eta;
p.Friction.Lhat_VW = Lhat_VW;
p.Friction.H_star_mean = H_s_M;
p.Friction.Symmetry = symmetry;
end
end