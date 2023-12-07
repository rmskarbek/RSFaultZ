function p = RSFaultZSetUp
%%% Get parameter values from fields. This function is called by
%%% the "Set Properties" button.
    GeometryOption = 'Thrust Fault';
    StateOption = 'Aging';
%%% Run time of the simulation [years].
    time = [0 100];
    
%%% Frictional parameters.
    a = 0.01;
    ab = 1.25e-4;
    b = a - ab;
    d_c = 2e-4;
    v_ref = 1e-9;
    mu_ref = 0.3;

    v_plate = 1e-9;

%%% Elastic parameters
    v_shear = 3000;                         % shear wave speed [m / s]
    Poisson = 0.25;                         % Poisson ratio
    ShearMod = 3e3;
    eta = ShearMod/(2*v_shear);             % radiation damping [MPa s / m]
    EffStress_0 = 10;                       % background normal stress [MPa]
    
    
    switch(GeometryOption)
        case 'Full Space'
            Poisson_1 = app.PoissonRatioEditField.Value;
            ShearMod_1 = 1e3*app.ShearModulusGPaEditField.Value;
            ShearMod_2 = ShearMod_1;
            Poisson_2 = Poisson_1;
            
            Beta = -(ShearMod_2*(1 - 2*Poisson_1) - ShearMod_1*(1 - 2*Poisson_2))...
                ./(2*(ShearMod_2*(1 - Poisson_1) + ShearMod_1*(1 - Poisson_2)));
            
            M = (2*ShearMod_1*ShearMod_2/(1 - Beta^2))/((ShearMod_2*(1 - Poisson_1)...
                + ShearMod_1*(1 - Poisson_2)));
            eta = ShearMod_1/(2*v_shear);             % radiation damping [MPa s / m]
            
%%% Length scales.    
            H_s = pi*M*d_c/(EffStress_0*(b-a)); 
            L_b = d_c*M/(b*EffStress_0);
            L_c = 4*pi*d_c*M/(EffStress_0*a*mu_ref^2);
        
%%% Set up the numerical grid.
            FaultLength = 1e3*app.FaultLengthkmEditField.Value;
            dxi = L_b/40;
            N_40 = round(FaultLength/dxi);
            n2 = nextpow2(N_40);
            if 2^n2 < 0.9*N_40
                N = 2^(n2 + 1);
            else
                N = 2^n2;
            end
            N_tot = 3*N;
            dx = FaultLength/(N - 1);
            TotalLength = N_tot*dx;
            %TotalLength = FaultLength + 2*PadLength;
            Xi = (0 : dx : FaultLength)' - FaultLength/2;
            X_tot = (0 : dx : TotalLength)' - floor(N_tot/(2*N))*N*dx;
            N_pad = (N_tot - N)/2;
            
%%% Computate the wave number.
            WaveNumber = fftshift(2*pi.*abs((-N_tot/2:N_tot/2-1)')./(TotalLength-dx));
            
        case {'Thrust Fault', 'Normal Fault'}                     
%%% Length scales.    
            H_s = pi*ShearMod*d_c/(EffStress_0*(1 - Poisson)*(b-a));
            L_b = d_c*ShearMod/(b*EffStress_0*(1 - Poisson));
            L_c = 4*pi*d_c*ShearMod/(EffStress_0*(1 - Poisson)*a*mu_ref^2);
            
%%% Set up grid: xi is the along-fault coordinate, x is the horizontal
%%% coordinate, z is the vertical coordinate. Grid is defined on segment 
%%% centers.
            FaultLength = 50e3;

            % dxi = 250;
            % N = round(FaultLength/dxi);

            N = 2^11;
            dxi = FaultLength/(N - 1);

            Xi = 0*dxi + (dxi/2: dxi: N*dxi - dxi/2)';

%%% Wedge geometry.
            alpha = 0;
            beta = 5;
%%% Zed is the fault's depth below the upper surface.
            Zed = Xi*(sin(pi*beta/180) + cos(pi*beta/180)*tan(pi*alpha/180));
%%% Horizontal location.
            X = Xi*cos(pi*beta/180);               % fault depth [m]
            
%%% Compute the elastic kernels for the defined grid.
            K_s = hodlr('handle',...
              @(i,j) WedgeKernel1b_HODLR(i, j, Xi', beta, alpha, 1), N, N); 
            K_n = hodlr('handle',...
             @(i,j) WedgeKernel1b_HODLR(i, j, Xi', beta, alpha, 0), N, N);
%                 K_s = full(K_s);
%                 K_n = full(K_n);
    end

%%% Loading rate.
    Loading = 0;            
    % Loading = (v_plate/pi)*(pi/2 ...
    %     + asin((FaultLength/2 - (Xi - Xi(1)))/(FaultLength/2))) +1e-3*v_plate;
    
%%% Initial slip velocity.
    v_0 = 1e0*v_plate*ones(size(Xi));
%%% Initial normal stress.
    sigma_0 = EffStress_0.*ones(size(Xi));
%%% Initial shear stress.
    tau_0 = mu_ref*sigma_0;
%%% Initial state.
    state_0 = 9e-1*(d_c./v_plate)*ones(size(Xi));

%%% Add load
    %v_0 = v_0 + v_load;
    Vars_0 = [v_0; tau_0; state_0; sigma_0];
    
%%% Friction parameters.
    p.Friction = struct('a', a, 'b', b, 'd_c', d_c, 'v_ref', v_ref,...
        'mu_ref', mu_ref, 'H_star', H_s, 'L_b', L_b, 'L_c', L_c);
%%% Options.
    plots = true;
    p.Options = struct('StateLaw', StateOption, 'Geometry', GeometryOption,...
        'Restart', 0, 'Plots', plots, 'Import', 1);

%%% Different parameters for the different geometries.
    switch(GeometryOption)
        case 'Full Space'
%%% Geometrical and numerical parameters.    
            p.Geometry = struct('FaultLength', FaultLength,...
                'FaultCoordinates', Xi, 'GridPoints', N,...
                'GridSpacing', dxi, 'Padding', N_pad, 'WaveNumber', WaveNumber,...
                'Runtime', time, 'Loading', Loading, 'InitConditions', Vars_0);
%%% Material parameters.
            p.Material = struct('ElasticModulus', M, 'BiMatContrast', Beta,...
                'BackgroundStress', EffStress_0, 'PlateRate', v_plate, 'ShearWaveSpeed',...
                v_shear, 'RadiationDamping', eta);
            
        case {'Thrust Fault', 'Normal Fault'}            
%%% Geometrical and numerical parameters.    
            p.Geometry = struct('SurfaceSlope', alpha, 'FaultDip', beta, 'FaultLength',...
                FaultLength, 'FaultCoordinates', Xi, 'GridPoints', N, 'GridSpacing',...
                dxi, 'Runtime', time, 'Loading', Loading, 'InitConditions', Vars_0,...
                'ShearKernel', K_s, 'NormalKernel', K_n, 'Depth', Zed,...
                'Distance', X);            
%%% Material parameters.
            p.Material = struct('ShearModulus', ShearMod, 'PoissonRatio', Poisson,...
                'BackgroundStress', EffStress_0, 'PlateRate', v_plate, 'ShearWaveSpeed',...
                v_shear, 'RadiationDamping', eta);
    end
    
end