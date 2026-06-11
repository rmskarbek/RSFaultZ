function Kernel = WedgeStress_HODLR_LS3(i, j, Z, Z_0, alpha, flag, type, sense, varargin)

%%% 3.5.2026 - The extra input is Wp - the complex derivative of the parameterized
%%% coordinates, at the cell edges. It is unused, so I'm not sure of the intended
%%% purpose. Probably, I stopped working on this problem before I implemented
%%% whatever I was thinking of.

%%% This code follows the notation in Skarbek (2025).
% Z     - Complex coordinates of recievers, e.g. fault centers, or anywhere else
% Z_0   - Complex coordinates of dislocations (i.e. sources), e.g. fault cell edges.

%%% This code assumes that slip and slip rate are positive. Because of this, the sign
%%% of the normal stress kernel is changed depending on the sense of slip (thrust, or
%%% normal fault) and the dip angle (beta) of the fault. Important: the dip angle is 
%%% measured from the x-axis towards the y-axis.

%%% For the normal stress kernel:
%%% 1. Thrust Fault, beta < 90: Kernel --> -Kernel 

%%% 2. Thrust Fault, beta > 90: no change             

%%% 3. Normal Fault, beta < 90: no change

%%% 4. Normal Fault, beta > 90: Kernel --> -Kernel


%%%------------------------------------------------------------------------------------%%%
%%% Compute the local tangent angles at the cell centers and edges.
%%%------------------------------------------------------------------------------------%%%
%%% Compute local tangent (or dip) angles at the cell centers. These angles will be used
%%% to rotate the coordinates for computing shear and normal stresses at the fault 
%%% centers. If stresses were desired at locations other than on the fault, then Theta
%%% would need to be altered accordingly.
    dx = real(Z_0(2:end) - Z_0(1:end-1));
    dy = imag(Z_0(2:end) - Z_0(1:end-1));

%%% 7.18.2025 - The absolute value will screw things up if the fault dip is greater
%%% than pi.
    % Theta = abs(atan(dy./dx));                                              % [rad]
    Theta = atan(dy./dx);                                                   % [rad]
    % R = pi;
    % Theta = Theta + R;

%%% Compute the local tangent angle at the cell edges. These angles determine the
%%% direction of the slip vector on the fault. Here we assume that slip only occurs
%%% parallel to the fault, i.e. there is no opening/closing component of slip.

%%% This implements the case illustrared by "Discretized fault for full elastic response"
%%% in Figure F1 of Romanet et al., GJI (2020). Each cell's edges are assigned the same
%%% tangent angle that occurs at that cell center. This means that edges connecting
%%% adjacent cells have two tangent angles, one for each cell on either side of the edge.
    Beta = Theta;                                % [rad]

%%% This implements the case illustrared by "Discretized fault curvature term only" in
%%% Figure F1 of Romanet et al., GJI (2020).
    switch type
        case 'single'
            Beta = [Beta(1); movmean(Beta, 2, 'EndPoints', 'discard'); Beta(end)];
            % Beta = Beta + R;
    end

%%% For accretionary wedge geometry, rotate beta and theta by the wedge surface slope.
    Beta = (pi/180)*alpha + Beta;
    Theta = (pi/180)*alpha + Theta;


%%%------------------------------------------------------------------------------------%%%
%%% Dislocation strength.
%%%------------------------------------------------------------------------------------%%%
%%% This assumes that slip is along the fault, in the direction of the dip angle. In other
%%% words, there is no component of slip that is normal to the fault surface. Gamma must 
%%% be transposed to correctly compute the complex potentials.
    Gamma = (-1i*exp(1i*Beta)/2).';
    

%%%------------------------------------------------------------------------------------%%%
%%% Additional transposes and indexing for correctly computing the complex potentials.
%%%------------------------------------------------------------------------------------%%%
    z = Z(i);
    Z_0 = Z_0.';
    theta = Theta(i);


%%%------------------------------------------------------------------------------------%%%
%%% For curved faults.
%%%------------------------------------------------------------------------------------%%%
%%% Provide the complex derivative of the parameterized coordinates, at the cell edges.
    if ~isempty(varargin)
        Wp = varargin{1,1};
    else
        Wp = ones(size(Beta));
    end


%%%------------------------------------------------------------------------------------%%%
%%% Evalute the potentials.
%%%------------------------------------------------------------------------------------%%%
    switch type
        case 'dipole'
%%% Up-dip edge of fault element.
            z_0 = Z_0(j);
            gamma = Gamma(j);
            wp = Wp(j);
            [Shear_1, Normal_1] = HalfSpacePotentials(z_0, z, theta, gamma, wp);
    
%%% Down-dip edge of fault element.
            z_0 = Z_0(j+1);
            gamma = Gamma(j);
            wp = Wp(j);
            [Shear_2, Normal_2] = HalfSpacePotentials(z_0, z, theta, gamma, wp);

        case 'single'
            % z = varargin{1,1};
            % theta = mean(angle(z));
            z_0 = Z_0(j);
            gamma = Gamma(j);
            wp = Wp(j);
            [Shear_1, Normal_1] = HalfSpacePotentials(z_0, z, theta, gamma, wp);
            % [Shear_2, Normal_2] = HalfSpacePotentials(z_0, z, theta, gamma);
        %%% 6.23.2025 - I THINK YOU INTENDED FOR Shear_1 and Normal_1 TO BE ZERO.
            Shear_2 = 0;
            Normal_2 = 0;
            % Shear_1 = 0;
            % Normal_1 = 0;
    end

%%% Assemble.
    if flag == 1
        Kernel = Shear_2 - Shear_1;
    else
        Kernel = Normal_2 - Normal_1;

%%% Deal with sign changes, as described above.
%%% NOTE 7.24.2025: Beta is a vector now, so there are issues here when the values of Beta
%%% are not all the same.
        switch sense
            case 'Thrust Fault'
                if Beta < pi/2
                    Kernel = -Kernel;
                end

            case 'Normal Fault'
                if Beta > pi/2
                    Kernel = -Kernel;
                end
        end
    end


%%%------------------------------------------------------------------------------------%%%
%%% Complex potentials for a dislocations in a half space. For example, see Skarbek, GJI
%%% (2025), Appendix C.
%%%------------------------------------------------------------------------------------%%%
    function [Shear, Normal] = HalfSpacePotentials(z_0, z, theta, gamma, wp)
%%% Transposed: z_0, gamma

%%% z_0 is the location of the dislocation. I.E. on-fault coordinate.
        zbar_0 = conj(z_0);
        gammabar = conj(gamma);

%%% Potentials. From equations (C13) and (C14) in Skarbek (2025).
%%% First derivative of big Omega.
        O_p = gamma.*(1./(z - z_0) - 1./(z - zbar_0))...
            - gammabar.*(z_0 - zbar_0)./(z - zbar_0).^2;

%%% Second derivative of big Omega.
        O_pp = gamma.*(1./(z - zbar_0).^2  - 1./(z - z_0).^2)...
            + 2*gammabar.*(z_0 - zbar_0)./(z - zbar_0).^3;
        
%%% First derivative of small omega.
        w_p = gammabar.*(1./(z - z_0) - 1./(z - zbar_0))...
            + gamma.*zbar_0.*(1./(z - z_0).^2 - 1./(z - zbar_0).^2)...
            + gammabar.*(zbar_0 - z_0).*(z + zbar_0)./(z - zbar_0).^3;        

%%% Complex Traction.
        % Traction = wp.'.*(O_p + conj(O_p) + exp(-2i*theta).*(z.*conj(O_pp) + conj(w_p)));
        Traction = (O_p + conj(O_p) + exp(-2i*theta).*(z.*conj(O_pp) + conj(w_p)));

%%% Shear stress due to shear dislocation.
        Shear = -imag(Traction);
        
%%% Normal stress due to shear dislocation.
        Normal = real(Traction);
    end

end