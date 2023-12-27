function BasicMovie(Out)
%%% This function provides some basic commands to animate the results of a
%%% simulation.

%Out = SimData2;
% FaultLength = Out.Params.Geometry.FaultLength;
% N = Out.Params.Geometry.GridPoints;
% N_t = numel(Out.Output.Time);
% beta = Out.Params.Geometry.FaultDip;
% Z = Out.Params.Geometry.Depth;

FaultLength = Out.Params.Geometry.FaultLength;
N = Out.Params.Geometry.GridPoints;
N_t = numel(Out.Output.Time);

X = Out.Params.Geometry.FaultCoordinates;
% beta = Out.Params.Geometry.FaultDip;
% Z = Out.Params.Geometry.Depth;
% X = Z./tan(pi*beta/180);

frames = floor(N_t/10);
F(frames) = struct('cdata',[],'colormap',[]);

Vel = Out.Output.Vel;


t = figure;
set(gcf,'Color','w')
i_s = 85;
spy = 365.25*24*3600;
s = 50;
for i = 1:frames
    semilogy(X, Vel(i*s, :), 'b')
    time = Out.Output.Time(i*s);
    
    r = mod(time(1), spy);
    years = (time(1) - r)/spy;
    days = r/(24*3600);
    caption = sprintf('%8d years, %6.2f days', [years, days]);
    title(caption);
        
    % xlim([0 round(max(X))])
    ylim(10.^[-13 0])
    xlabel('Distance (m)')
    ylabel('Slip Velocity (m / s)')
%    set(gca, 'Xdir', 'reverse')
    drawnow
    F(i) = getframe(gcf);
end

% v = VideoWriter('movie.avi');
% open(v);
% writeVideo(v,F);
% close(v)