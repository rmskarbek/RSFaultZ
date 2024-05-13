function SourceFunction = SinusoidSource1(Xi)
%%% Time-dependent sinusoid fluid pressure source function. Does not depend on along-fault
%%% distance.

%%%--------------------------------------------------------------------------------------%%%
%%%--------------------------------------------------------------------------------------%%%
%%% Flow flow parameters from BP6.
    PeriodSource = 24*3600*100;          % [s]
    AmplitudeSource = 10;               % [MPa]

%%% Apply the same amplitude to the entire fault.
    AmplitudeSource = AmplitudeSource*ones(numel(Xi),1);
    
%%% Create the source function.   
    SourceFunction = @SourceFunc;

    function SourceTerm = SourceFunc(t)

        SourceTerm = (2*pi*AmplitudeSource/PeriodSource)*cos(2*pi*t/PeriodSource);

    end

end
