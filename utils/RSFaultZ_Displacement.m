function [UX, UY, Time, K_x, K_y] = RSFaultZ_Displacement(SimData, C, T)

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%%% This function computes displacements from RSFaultZ simulation output.

%%% Inputs
% SimData   - structure output file from RSFaultZ

% C         - [2 x M] array of coordinates of M locations where displacement
%             is desired. Column 1: x coordinates. Column 2: y coordinates.
%             Units [km]

% T         - Optional downsampling index for output times. Default value
%             is 1. E.G. if T = 10, displacements are computed for every
%             tenth time in the the simulations output.

%%% Outputs
% UX        - Array of x-displacements. Each row gives the displacement for
%             cooresponding row of C.

% UY        - Same as UX, for y-displacements.

%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

%%% Fault coordinates.
Xi = SimData.Params.Geometry.FaultCoordinates;
N = SimData.Params.Geometry.GridPoints;

%%% Geometry.
geom = SimData.Params.Options.Geometry;
beta = SimData.Params.Geometry.FaultDip;
alpha = SimData.Params.Geometry.SurfaceSlope;

%%% Elastic constants.
ShearMod = SimData.Params.Material.ShearModulus;
Poisson = SimData.Params.Material.PoissonRatio;

%%% Convert displacement locations to complex coordinates.
Z = complex(C(:,1), C(:,2));
M = numel(Z);

%%% Compute the elastic kernels for the list of coordinates. The kernels 
%%% will have dimensions [M x N], where N is the number of grid points on
%%% the fault.

%%% For now, form the full matrices by evaluating one row at a time.
K_x = nan(M, N);
K_y = nan(M, N);
j = 1:N;
for q = 1:M
    Kx = WedgeDisp_HODLR(j, q, Xi', Z', beta, alpha, ShearMod, Poisson, 1);
    Ky = WedgeDisp_HODLR(j, q, Xi', Z', beta, alpha, ShearMod, Poisson, 0);
    K_x(q,:) = Kx;
    K_y(q,:) = Ky;
end

%%% Get times and slips for computing displacements.
Time = downsample(SimData.Output.Time, T);
N_t = numel(Time);
Slip = downsample(SimData.Output.Slip, T)';

%%% Allocate displacement arrays
UX = nan(M, N_t);
UY = nan(M, N_t);

%%% Compute displacements. For normal faults the sense of slip is opposite
%%% from thrust faults, so the sign of all displacements changes.
for k = 1:N_t
    switch geom
        case 'Thrust Fault'
            UX(:,k) = K_x*Slip(:,k);
            UY(:,k) = K_y*Slip(:,k);
        case 'Normal Fault'
            UX(:,k) = -K_x*Slip(:,k);
            UY(:,k) = -K_y*Slip(:,k);
    end
end