function [T_n, D_n, V_m, V_max, Time] = FaultZAnalyze1(SimData)

spy = 365.25*24*3600;
Time = SimData.Output.Time/spy;
[V_m, I] = max(SimData.Output.Vel,[],2);
dV_m = diff(V_m);
N_t = numel(V_m);
V_max = max(V_m(round(N_t/2):N_t));

d = V_m > 1e-3;
if isempty(V_m(d)) == 0
    k = 1;
    for i = 1:numel(dV_m)
        if V_m(i) < 1e-3 && V_m(i) + dV_m(i) > 1e-3
            if abs(V_m(i) - 1e-3) < abs(V_m(i+1) - 1e-3)
                j(k) = i; 
                k = k + 1;
            else
                j(k) = i+1; 
                k = k+1;
            end
        end
    end
    D_n = SimData.Params.Geometry.Distance(I(j))/1e3;
    T_n = SimData.Output.Time(j)/spy;
else
    D_n = [];
    T_n = [];
end