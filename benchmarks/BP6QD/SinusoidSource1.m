function SourceFunction = SinusoidSource1(Xi)
%%% Time-dependent sinusoid fluid pressure source function. Does not depend on along-fault
%%% distance.

%%%--------------------------------------------------------------------------------------%%%
%%%--------------------------------------------------------------------------------------%%%
%%% Flow flow parameters from BP6.
    % spd = 24*3600;
    % PeriodDays = 365.25;    
    % PeriodSource = spd*PeriodDays;        % [s]

    spy = 365.25*24*3600;
    PeriodYears = 1;
    PeriodSource = spy*PeriodYears;        % [s]
    AmplitudeSource = 0.0217/2;             % [MPa]

%%% Apply the same amplitude to the entire fault.
    AmplitudeSource = AmplitudeSource*ones(numel(Xi),1);
    
%%% Create the source function.   
    SourceFunction = @SourceFunc;

    function SourceTerm = SourceFunc(t)

        SourceTerm = (2*pi*AmplitudeSource/PeriodSource)*cos(2*pi*t/PeriodSource);

    end

end
