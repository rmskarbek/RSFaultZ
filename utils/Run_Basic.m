%%% This is a basic script to programatically run RSFaultZ.

%%% Create an instance of the app.
R = RSFaultZ;

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Set the default properties.
R.SetPropertiesButton.ButtonPushedFcn([],[]);

%%% Change some options.
R.StateLawDropDown.Value = 'Aging';
% R.PlotsOnCheckBox.Value = false;

%%% Run the simulation.
R.StartButton.ButtonPushedFcn([],[]);

%%% Access the simulation output.
SimData = struct('Params', R.p, 'Output', R.Out);

%%% Delete the app instance.
delete(R);