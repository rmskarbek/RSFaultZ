function Kernel = WedgeStress_HODLR(i, j, xi, beta, alpha, d, flag, varargin)

%%% This file needs to be cleaned up. I'm not sure if the observation vector option is still working.
  
%%% To evaluate for full matrix, where x = xi are Nx1 vectors:
%%% Kernel = WedgeStress_HODLR((1:N), (1:N), Xi', beta, alpha, 1);

%%% To evaluate for full matrix, where x is Mx1, and xi is Nx1, and no
%%% element of z is equal to any element of xi:
%%% Kernel = WedgeStress_HODLR((1:M), (1:N), Xi', beta, alpha, 1);

%%% Direction of the burger's vector, i.e. slip is along a fault dippping
%%% at an angle beta.
    beta = pi - (alpha + beta)*pi/180;
%%% Dislocation strength, pointing along the fault at angle beta.
    gamma = -1i*exp(1i*beta)/2;
    

%%% Locations where stress is desired are z. Dislocations are located at s.    
    if isempty(varargin)
%%% If no observation vector is given, stresses are computed on the fault.      
%%% Theta is the angular coordinate of observation points, which are on the
%%% fault.
        theta = beta;
        % z = xi(i)'*exp(1i*theta);
        z = complex(xi(i)'*cos(beta), d + xi(i)'*sin(beta));

%%% Grid spacing.
        dxi = mean(diff(xi));

%%% Up-dip edge of fault element.
        s = xi(j) - dxi/2;
        [Shear_1, Normal_1] = HalfSpacePotentials(s, z, beta, theta, gamma);
    
%%% Down-dip edge of fault element.
        s = xi(j) + dxi/2;
        [Shear_2, Normal_2] = HalfSpacePotentials(s, z, beta, theta, gamma);
    else
%%% If an observation vector is supplied, stresses are comuputed at those
%%% coordinates.
        z = varargin{1,1};
        theta = mean(angle(z));
        s = xi(j);
        [Shear_1, Normal_1] = WedgePotentials(s, z, q, beta, theta, gamma);
        Shear_2 = 0;
        Normal_2 = 0;
    end

%%% Assemble.
    if flag == 1
        Kernel = Shear_2 - Shear_1;
    else
        Kernel = Normal_2 - Normal_1;
    end               


    function [Shear, Normal] = HalfSpacePotentials(s, z, beta, theta, gamma)
        % gamma = -(1i/2)*complex(cos(beta), sin(beta));
        gammabar = conj(gamma);

%%% z_0 is the location of the dislocation. I.E. on-fault coordinate.
        % z_0 = s*exp(1i*beta);
        z_0 = complex(s*cos(beta), d + s*sin(beta));
        zbar_0 = conj(z_0);

%%% Potentials. BF1, p. 54b; BF2 p.90b
        F = gamma.*(1./(z - z_0) - 1./(z - zbar_0))...
            - gammabar.*(z_0 - zbar_0)./(z - zbar_0).^2;

        Fp = gamma.*(1./(z - zbar_0).^2  - 1./(z - z_0).^2)...
            + 2*gammabar.*(z_0 - zbar_0)./(z - zbar_0).^3;
        
        Y = gammabar.*(1./(z - z_0) - 1./(z - zbar_0))...
            + gamma.*zbar_0.*(1./(z - z_0).^2 - 1./(z - zbar_0).^2)...
            + gammabar.*(zbar_0 - z_0).*(z + zbar_0)./(z - zbar_0).^3;        

%%% Complex Traction.
        Traction = F + conj(F) + exp(-2i*theta).*(z.*conj(Fp) + conj(Y));
%%% Shear stress due to shear dislocation
        Shear = -imag(Traction);
%%% Normal stress due to shear dislocation.
        Normal = real(Traction);
    end

end
