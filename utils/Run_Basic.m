%%% This is a basic script to programatically run RSFaultZ.

%%% Create an instance of the app.
R = faultZ;

%%% Set the default properties.
R.SetPropertiesButton.ButtonPushedFcn([],[]);

%%% Change some options.
R.StateLawDropDown.Value = 'Aging';
R.PlotsOnCheckBox.Value = false;

%%% Plot the geometry.
R.PlotGeometryButton.ButtonPushedFcn([],[]);

%%% Run the simulation.
R.StartButton.ButtonPushedFcn([],[]);

%%% Acess the simulation output.
Out = R.Out;