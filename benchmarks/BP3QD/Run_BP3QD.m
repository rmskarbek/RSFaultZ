function SimData = Run_BP3QD(geom, DipAngle)

%%% This script will programatically run the SCEC benchmark problem BP3QD.

%%% Execute time ~4.5 minutes for 60 deg thrust fault.
%%% Execute time ~3.6 minutes for 30 deg normal fault.

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%% Geometry.
FaultLength = 40;

%%% Friction parameters.
a = 0.01;
b = 0.015;
d_c = 0.008;
v_ref = 1e-6;
mu_ref = 0.6;

%%% Elastic parameters.
v_shear = 3464;
density = 2670;
Poisson = 0.25;
ShearMod = v_shear^2*density/1e9;

%%% Background effective stress.
EffStress_0 = 50;

%%% Plate rate.
v_plate = 1e-9;

%%% Run time.
RunTime = 700;

%%% Grid spacing.
dxi = 25;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%%% Create an instance of the app.
R = RSFaultZ;

%%% Set options.
R.StateLawDropDown.Value = 'Aging';
R.PlotsOnCheckBox.Value = false;
    
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

%%% Set run time.
R.RunTimeyrEditField.Value = RunTime;
    
%%% Set nominal frictional parameters.
R.aEditField.Value = a;
R.bEditField.Value = b;
R.dcmEditField.Value = d_c;
R.v0EditField.Value = v_ref;
R.mu0EditField.Value = mu_ref;

%%% Set plate rate.
R.LoadingRatemsEditField.Value = v_plate;

%%% Set grid spacing.
R.dximEditField.Value = dxi;

%%% Update the Grid Controls panel.
R.dximEditField.ValueChangedFcn([],[]);

%%% Set the properties to compute grid coordinates.
R.SetPropertiesButton.ButtonPushedFcn([],[]);

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Get fault coordinates.
N = R.p.Geometry.GridPoints;
Xi = R.p.Geometry.FaultCoordinates;

%%% Find coordinates of transitions in RSF parameter a.
[~, i2] = min(abs(Xi - 18e3));
[~, i1] = min(abs(Xi - 15e3));

%%% Compute depth-varying a.
a_max = 0.025;
a = a*ones(N, 1);
a(i1:i2) = linspace(0.01, a_max, numel(a(i1:i2)));
a(i2:end) = a_max;

%%% Set depth-varying a.
R.p.Friction.a = a;

%%% Get radiation damping term.
eta = R.p.Material.RadiationDamping;

%%% Initial slip velocity.
v_i = 1e0*v_plate*ones(size(Xi));
%%% Initial normal stress.
sigma_i = EffStress_0.*ones(size(Xi));
%%% Initial shear stress.
tau_i = EffStress_0.*a_max*asinh((v_i./(2*v_ref)).*exp((mu_ref...
    + b*log(v_ref./v_i))./a_max)) + eta*v_i;
%%% Initial state.
state_i = (d_c/v_ref)*exp((a./b).*log((2*v_ref./v_i)...
    .*sinh((tau_i - eta*v_i)./(a*EffStress_0))) - mu_ref/b);

%%% Set initial conditions.
Vars_i = [v_i; tau_i; state_i; sigma_i];
R.p.Geometry.InitConditions = Vars_i;


%%% Set import flag to 1, to prevent altered parameters from being
%%% overwritten.
R.p.Options.Import = 1;

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Run the simulation.
R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
SimData = struct('Params', R.p, 'Output', R.Out);

%%% Close the RSFaultZ instance.
delete(R);
end