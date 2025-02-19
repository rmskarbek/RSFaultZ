function [R, p] = StabilitySetup(ab, Length_hat, beta, mu_ref, geom,...
    Burial_hat, type, varargin)

%%% This function creates an instance of RSFaultZ, or alters an already
%%% existing instance, and sets parameters for numerical stability
%%% analysis. It is called by StabilityRun1.m

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%% Friction parameters.
b = 0.01;
%%% Value of a for VW.
a = b*ab;
%%% Value of a for VS.
% a = ab + b;

%%% d_c for linear stability calculations.
d_c = 0.002;
%%% RSF reference velocity.
v_ref = 1e-6;

%%% Elastic parameters.
v_shear = 3000;
Poisson = 0.25;
%%% Shear modulus in GPa.
ShearMod = 30;
M = 1e3*ShearMod/(1 - Poisson);

%%% Background effective stress.
EffStress_0 = 5;

%%% Friction length scales.
H_s = pi*M*d_c/(EffStress_0*(b-a));
L_b = d_c*M/(b*EffStress_0);

%%% For thin layer geometry, the length scales are modified.
switch geom
    case 'Layer'
        d = Burial_hat*L_b;            % [m]
        H_s = (8*pi*d*H_s).^(1/2);
        L_b = (2*d*L_b).^(1/2);
        BurialDepth = d/1e3;            % [km]
end

%%% Set grid fault length, grid spacing, and burial depth.
if a - b < 0

%%% Fault length [km] and grid spacing [m] for VW.
    FaultLength = Length_hat*H_s/1e3;

    switch geom
        case {'Thrust Fault', 'Normal Fault'}
    %%% For dipping critical fault length, make sure there are at least 250
    %%% grid points.        
            dxi1 = L_b/80;
            dxi2 = 1e3*FaultLength/250;
            if dxi2 < dxi1
                dxi = dxi2;
            else
                dxi = dxi1;
            end
            BurialDepth = Burial_hat*H_s/1e3;  % [km]

        case {'Full Space', 'Layer'}
%%% For full-space and thin layer calculations
            switch type
                case 'search'
                    dxi = L_b/80;
                case 'wavelength'
                    dxi = L_b/20;
            end
    end      
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
    
%%% Set geometry. Use the 'Full Space' geometry for thin layer because it's faster. The
%%% stress change functions for the thin layer geometry will be defined in
%%% StabilityAnalysis.m.
if strcmp(geom, 'Layer') == 1
    R.GeometryDropDown.Value = 'Thrust Fault';
else
    R.GeometryDropDown.Value = geom;
end

R.FaultLengthkmEditField.Value = FaultLength;
if strcmp(geom, 'Full Space') == 0
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

%%% Set run time.
R.RunTimeyrEditField.Value = 80;

%%% Set the properties to compute grid coordinates.
R.SetPropertiesButton.ButtonPushedFcn([],[]);

%%% Set import flag to true, to prevent altered parameters from being
%%% overwritten.
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