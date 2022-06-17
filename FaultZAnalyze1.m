function [T_n, D_n] = FaultZAnalyze1(SimData)

spy = 365.25*24*3600;
[V_m, I] = max(SimData.Output.Vel,[],2);
dV_m = diff(V_m);

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