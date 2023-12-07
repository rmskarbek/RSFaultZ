function [R, p] = RSFaultZ_setup

%%% This script will programatically set up a basic simulation.

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%% Geometry.
FaultLength = 40;
DipAngle = 60;

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
R = faultZ;

%%% Set options.
R.GeometryDropDown.Value = 'Full Space';
R.StateLawDropDown.Value = 'Aging';
R.PlotsOnCheckBox.Value = false;
    
%%% Set geometry.
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
    
%%% Set frictional parameters.
R.aEditField.Value = a;
R.bEditField.Value = b;
R.dcEditField.Value = d_c;
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

%%% NEED TO SET import flag TO TRUE OR WHATEVER, TO PREVENT ALTERED VECTORS
%%% FROM BEING OVERWRITTEN.
R.p.Options.Import = 1;

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Get parameters structure.
p = R.p;

%%% Run the simulation.
% R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
% SimData = struct('p', R.p, 'Output', R.Out);

end