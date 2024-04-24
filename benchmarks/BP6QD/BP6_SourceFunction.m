function SourceFunction = BP6_SourceFunction(Xi)
%%% This function computes the source term for benchmark problem BP6QD. To evaluate the
%%% dirac delta function, it makes use of:

%%% Divyaprakash (2024). Smoothed Dirac Delta Function.
%%% (https://www.mathworks.com/matlabcentral/fileexchange/157136-smoothed-dirac-delta-function),
%%%% MATLAB Central File Exchange. Retrieved April 8, 2024.

%%%--------------------------------------------------------------------------------------%%%
%%%--------------------------------------------------------------------------------------%%%
%%% Flow flow parameters from BP6.
    TimeSource = 24*3600*100;           % [s]
    Compress = 1e6*1e-8;                % [1 / MPa]
    phi_0 = 0.1;                            
    FlowRate = 1.25e-6;                  % [m / s]

%%% Numerically evaluate the Dirac delta function for the simulation grid.    
    dxi = mean(diff(Xi));    
    Dirac = EvaluateDirac(dxi, Xi);
    FlowRate = FlowRate*Dirac;

%%% Pressure source.
    Source = FlowRate/(phi_0*Compress);
    SourceFunction = @SourceFunc;

    function SourceTerm = SourceFunc(t)

        if t < TimeSource
            SourceTerm = Source;
        else
            SourceTerm = 0;
        end

    end

    function Z = EvaluateDirac(dx, X)
    % Evaluate the dirac delta function at a point in the middle of the domain
        xq = X(round(end/2),1);
        Z = zeros(size(X)); 
        
    % Loop through each point in the grid and evaluate the dirac delta function
        for i = 1:size(X, 1)
        % Compute the distance vector from the current point to the query point
            xin = X(i,1) - xq;
            
        % Evaluate the dirac delta function using the defined diracdelta function. Need to
        % rescale for 1D delta function.
            Z(i,1) = diracdelta(xin, dx)*dx;
        end
    end

end