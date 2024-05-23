function [C, Time, T_UX, T_UY] = BP3QD_Displacements(SimData)
%%% NOTES
%%% Rewrite this code into similar structure of OnFault code. Generate a
%%% structure with a field for each station location, that contains a table
%%% of the displacements. This structure can be used for writing text files
%%% that can be submitted as benchmark data.

%%% Off-fault station locations in complex polar coordinates.
Rho = 1e3*[(-32:8:-8)'; -1e-6; 1e-6; (8:8:32)'];
Theta = [pi*ones(numel(Rho)/2,1); zeros(numel(Rho)/2,1)];
Z = Rho.*exp(1i*Theta);

C = [Rho, zeros(numel(Rho),1)];
T = 10;

%%% The rigid body displacements and potentially also the back slips need
%%% to be redefined to account for the opposite sense of slip for normal
%%% and thrust faults. This is accounted for in RSFaultZ_Displacement.m,
%%% but not here. The correct displacements are computed here for
%%% thrust faults without any sign changes. 

%%% To obtain the correct displacements for normal faults, the sign of
%%% output from RSFaultZ_Displacement.m must be changed. Then the sign must
%%% be changed again after applying the back slips and rigid body
%%% displacements.
[UX, UY, Time, K_x, K_y] = RSFaultZ_Displacement(SimData, C, T);
geom = SimData.Params.Options.Geometry;
switch geom
    case 'Normal Fault'
            UX = -UX;
            UY = -UY;
end

%%% Plate rate.
N = SimData.Params.Geometry.GridPoints;
FaultDip = SimData.Params.Geometry.FaultDip;
v_plate = SimData.Params.InitialConditions.PlateRate;

%%% Rigid body displacements.
UX_fw = (v_plate*Time'/2)*cosd(FaultDip);
UX_hw = -(v_plate*Time'/2)*cosd(FaultDip);

UY_fw = (v_plate*Time'/2)*sind(FaultDip);
UY_hw = -(v_plate*Time'/2)*sind(FaultDip);


%%% Compute back slips.
BackSlip  = Time*v_plate;
UX_b = nan(size(C, 1), numel(Time));
UY_b = nan(size(C, 1), numel(Time));
for q = 1:numel(Time)
    UX_b(:,q) = K_x*BackSlip(q)*ones(N,1);
    UY_b(:,q) = K_y*BackSlip(q)*ones(N,1);
end

%%% Apply BP3QD far-field boundary conditions, assuming fault dip is
%%% 0 - 90 deg.
for q = 1:size(C,1)
    if C(q,1) < 0
        UX_rb = UX_hw;
        UY_rb = UY_hw;
    else
        UX_rb = UX_fw;
        UY_rb = UY_fw;
    end
    UX(q,:) = (UX_rb + UX_b(q,:)) - UX(q,:);
    UY(q,:) = (UY_rb - UY_b(q,:)) + UY(q,:);
end

switch geom
    case 'Normal Fault'
            UX = -UX;
            UY = -UY;
end

%%% Make a nice table. First create variable names for station locations.
%%% The x-axis in BP3-QD is defined in the opposite direction from that of
%%% RSFaultZ. The sign of the station locations is switched here to account
%%% for this.
VarNames = cell(size(C,1), 1);
for q = 1:size(C,1)
    if C(q,1) < 0
        VarNames{q,1} = ['p_',num2str(abs(C(q,1)/1e3)), 'km'];
    else
        VarNames{q,1} = ['n_',num2str(abs(C(q,1)/1e3)), 'km'];
    end
end
T_UX = array2table(UX', 'VariableNames', VarNames);
T_UY = array2table(UY', 'VariableNames', VarNames);