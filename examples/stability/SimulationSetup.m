function [R, p] = SimulationSetup(ab, Length_hat, beta, mu_ref, geom,...
    Burial_hat, EffStress_0, varargin)

%%% This function creates an instance of RSFaultZ, or alters an already
%%% existing instance, and sets parameters for full space simulations for
%%% confirmation of the critical fault length. It is called by
%%% SimulationRun1.m.

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%% Friction parameters.
b = 0.01;
%%% Value of a for VW.
a = b*ab;
%%% Value of a for VS.
% a = ab + b;

%%% d_c for full space simulations.
d_c = 0.01;
%%% RSF reference velocity.
v_ref = 1e-6;

%%% Elastic parameters.
v_shear = 3000;
Poisson = 0.25;
%%% Shear modulus in GPa.
ShearMod = 30;
M = 1e3*ShearMod/(1 - Poisson);

%%% Friction length scales.
H_s = pi*M*d_c/(EffStress_0*(b-a));
L_b = d_c*M/(b*EffStress_0);

%%% For thin layer geometry, the length scales are modified.
switch geom
    case 'Layer'
        d = Burial_hat*L_b;
        H_s = (8*pi*d*H_s).^(1/2);
        L_b = (2*d*L_b).^(1/2);
end

%%% Set grid fault length, grid spacing, and burial depth.
if a - b < 0

%%% Fault length [km] and grid spacing [m] for VW.
    FaultLength = Length_hat*H_s/1e3;

    switch geom
        case {'Thrust Fault', 'Normal Fault'}
%%% For dipping critical fault length, make sure there are at least 250
%%% grid points.
            dxi = 1e3*FaultLength/250;

        case {'Full Space', 'Layer'}
%%% For full-space and thin layer calculations.
            dxi = L_b/80;
            % dxi = 1e3*FaultLength/250;
    end

    %%% Units might be wrong here.
    BurialDepth = Burial_hat*H_s/1e3;
else

%%% Fault length [km] and grid spacing [m] for VS.
    FaultLength = Length_hat*L_b/1e3;
    dxi = 1e3*FaultLength/500;
    BurialDepth = 1e-3*M*d_c*Burial_hat/EffStress_0;
end
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

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
R.PlotsOnCheckBox.Value = false;
    
%%% Set geometry.
if strcmp(geom, 'Layer') == 1
    R.GeometryDropDown.Value = 'Full Space';
else
    R.GeometryDropDown.Value = geom;
end

R.FaultLengthkmEditField.Value = FaultLength;
if strcmp(geom, 'Thrust Fault') == 1 || strcmp(geom, 'Normal Fault')
    R.DipAngleEditField.Value = beta;
    R.BurialDepthkmEditField.Value = BurialDepth;
end

%%% Set normal stress.
R.NormalStressMPaEditField.Value = EffStress_0;

%%% Set elastic parameters.
R.vsmsEditField.Value = v_shear;
R.PoissonRatioEditField.Value = Poisson;
R.ShearModulusGPaEditField.Value = ShearMod;
    
%%% Set frictional parameters.
R.aEditField.Value = a;
R.bEditField.Value = b;
R.dcmEditField.Value = d_c;
R.v0EditField.Value = v_ref;
R.mu0EditField.Value = mu_ref;

%%% Set grid spacing.
R.dximEditField.Value = dxi;

%%% Update the Grid Controls panel.
R.dximEditField.ValueChangedFcn([],[]);

%%% Set the run time. Depending on parameters, different run times are required to achieve a 
%%% stable limit cycle.
if ab == 0.1
    R.RunTimeyrEditField.Value = 200;
    
elseif ab == 0.3
    if Length_hat == 0.3725
        R.RunTimeyrEditField.Value = 600;
    else
        R.RunTimeyrEditField.Value = 200;
    end

elseif ab == 0.5
    if Length_hat == 0.375
        R.RunTimeyrEditField.Value = 900;    
    else
        R.RunTimeyrEditField.Value = 300;
    end

elseif ab == 0.7
    if Length_hat == 0.3725
        R.RunTimeyrEditField.Value = 2500;
    else
        R.RunTimeyrEditField.Value = 1200;
    end
end

%%% Increase run time when the perturbation is on the edge of the fault.
R.RunTimeyrEditField.Value = 1.5*R.RunTimeyrEditField.Value;

%%% Set the properties to compute grid coordinates.
R.SetPropertiesButton.ButtonPushedFcn([],[]);

%%% Set import flag to true, to prevent altered parameters from being
%%% overwritten.
R.p.Options.Import = 1;

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Get parameters structure.
p = R.p;
switch geom
    case 'Layer'
        p.Geometry.BurialDepth = d;
end

%%% Run the simulation.
% R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
% SimData = struct('p', R.p, 'Output', R.Out);

end