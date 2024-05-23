function R = Run_BP6QD

%%% This script will programatically run the SCEC benchmark problem BP6QD.

%---------------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------%
%%% Define the parameters values for the problem.

%%% Geometry.
geom = 'Full Space';
DipAngle = 0;
FaultLength = 40;

%%% Friction parameters.
a = 0.007;
b = 0.005;
d_c = 0.004;                            % [m]
v_ref = 1e-6;                           % [m/s]
mu_ref = 0.6;

%%% Elastic parameters.
v_shear = 3464;                         % [m/s]
density = 2670;                         % [kg / m^3]
Poisson = 0;
ShearMod = v_shear^2*density/1e9;       % [GPa]

%%% Fluid flow parameters.
phi_0 = 0.1;
Viscosity = 1e-9;                       % [MPa s]
Compressibility = 1e-2;                        % [1 / MPa]
Permeability = 1e-13;                           % [m^2]

%%% Name of fluid source function.
SourceFunctionName = 'BP6_SourceFunction';

%%% Background effective stress.
EffStress_0 = 50;                       % [MPa]

%%% Plate rate.
v_plate = 0;                            % [m/s]

%%% Initial velocity.
v_init = 1e-12;                         % [m/s]

%%% Run time.
RunTime = 2;                            % [years]

%%% Grid spacing.
dxi = 25;                               % [m]

%---------------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------%
%%% Create an instance of RSFaultZ and set the parameters values and initial conditions

%%% Create an instance of the app.
R = RSFaultZ;

%%% Set options.
R.StateLawDropDown.Value = 'Aging';
R.RunTimeyrEditField.Value = RunTime;

%---------------------------------------------------------------------------------------------%
%%% 1. Set the geometry, normal stress, elastic, friction, and fluid flow parameters.

%%% Set geometry.
R.GeometryDropDown.Value = geom;
R.FaultLengthkmEditField.Value = FaultLength;
R.DipAngleEditField.Value = DipAngle;

%%% Set normal stress.
R.NormalStressMPaEditField.Value = EffStress_0;

%%% Set elastic parameters.
R.vsmsEditField.Value = v_shear;
R.PoissonRatioEditField.Value = Poisson;
R.ShearModulusGPaEditField.Value = ShearMod;

%%% Enable fluid flow.
R.DiffusionCheckBox.Value = 1;
R.SourceFunctionCheckBox.Value = 1;
R.DiffusionCheckBox.ValueChangedFcn([],[]);
R.SourceFunctionCheckBox.ValueChangedFcn([],[]);

%%% Set fluid flow parameters.
R.SourceFunctionName.Value = SourceFunctionName;
R.PorosityEditField.Value = phi_0;
R.FluidViscosityMPasEditField.Value = Viscosity;
R.Compressibility1MPaEditField.Value = Compressibility;
R.Permeabilitym2EditField.Value = Permeability;
    
%%% Set frictional parameters.
R.aEditField.Value = a;
R.bEditField.Value = b;
R.dcmEditField.Value = d_c;
R.v0EditField.Value = v_ref;
R.mu0EditField.Value = mu_ref;

%---------------------------------------------------------------------------------------------%
%%% 2. Create the numerical grid.

%%% Set grid spacing.
R.dximEditField.Value = dxi;

%%% Update the Grid Controls panel.
R.dximEditField.ValueChangedFcn([],[]);

%%% Create the grid.
R.CreateGridButton.ButtonPushedFcn([],[]);

%---------------------------------------------------------------------------------------------%
%%% 3. Set the initial conditions.

%%% Get fault coordinates.
N = R.p.Geometry.GridPoints;
Xi = R.p.Geometry.FaultCoordinates;

%%% Set plate rate.
R.LoadingRatemsEditField.Value = v_plate;

%%% Select steady state shear stress intial condition. This won't do anything since the
%%% initial conditions are overwritten below.
R.SteadyStateConditionButtonGroup.SelectedObject.Value = false;

%%% Initial slip velocity.
R.InitialVelocitymsEditField.Value = v_init;
R.RandomVariationEditField.Value = 0;
v_i = v_init*ones(size(Xi));

%%% Get radiation damping term.
eta = R.p.Material.RadiationDamping;

%%% Initial normal stress.
sigma_i = EffStress_0.*ones(size(Xi));
%%% Initial shear stress.
tau_0 = 29.2;
tau_i = tau_0 + eta*v_i;
%%% Initial state.
state_i = (d_c/v_ref)*exp((a./b).*log((2*v_ref./v_i)...
    .*sinh((tau_0 - eta*v_i)./(a*EffStress_0))) - mu_ref/b);
%%% Initial pore pressure.
pressure_i = zeros(N,1);

%%% Set initial conditions.
Vars_i = [v_i; tau_i; state_i; sigma_i; pressure_i];
R.p.InitialConditions.Vars_init = Vars_i;

%%% Set import flag to 1, to prevent altered parameters from being overwritten.
R.p.Options.Import = 1;

%---------------------------------------------------------------------------------------------%
%%% 5. Plotting options.

%%% Turn off plots.
% R.PlotsOnCheckBox.Value = false;

%%% Plot the pore pressure
R.PlotTypeDropDown.Value = 'Pore Pressure';
R.PlotTypeDropDown.ValueChangedFcn([],[]);
R.yminEditField.Value = -1;
R.ymaxEditField.Value = 10;
R.PlotButton.ButtonPushedFcn([],[]);

%---------------------------------------------------------------------------------------------%
%---------------------------------------------------------------------------------------------%
%%% Run the simulation.
% R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
% SimData = struct('Params', R.p, 'Output', R.Out);

%%% Close the RSFaultZ instance.
% delete(R);
end