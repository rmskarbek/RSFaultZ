function OnFaultData = BP3QD_OnFault(SimData)

%%% Off-fault station locations. Distance down-dip (m)
C = 1e3*[(0:2.5:20)'; [25; 30; 35]];

%%% On-fault coordinates.
Xi = SimData.Params.Geometry.FaultCoordinates;

%%% This loop creates a data table for each station and puts it into a
%%% structure.
for q = 1:numel(C)

%%% First generate a name for the station.
    StationName = ['fltst_dp',num2str(abs(C(q,1)/1e2))];

%%% Find the grid index of the station.
    [~, i] = min(abs(Xi - C(q)));

%%% Get the simulation data for that index, and put it in a table for this
%%% station.
    StationData = [SimData.Output.Time, SimData.Output.Slip(:,i),...
        log10(SimData.Output.Vel(:,i)), SimData.Output.ShearStress(:,i),...
        SimData.Output.NormalStress(:,i), log10(SimData.Output.State(:,i))];

    VarNames = [{'Time_s'}, {'Slip_m'}, {'SlipRate_log10ms'},...
        {'ShearStress_MPa'}, {'NormalStress_MPa'}, {'State_log10s'}];
    
    DataTable = array2table(StationData, 'VariableNames', VarNames);

%%% Put the data table into a structure field named with the station name.
    OnFaultData.(StationName) = DataTable;
end