function Kernel = StrikeSlipStress_HODLR(i, j, xi, beta)

%%% Along-fault spatial coordinate.
beta = (pi/180)*beta;
theta = beta - pi/2;
z = xi(i)'*exp(1i*beta);

%%% Variables s are the integrands
dxi = mean(diff(xi));

%%% Build discretized elastic stress operator assuming piecewise continuous
%%% slip distribution.
%%% Edge 1.
s = xi(j) - dxi/2;
Shear_1 = StrikeSlipPotentials(s, z, beta, theta);

%%% Edge 2.
s = xi(j) + dxi/2;
Shear_2 = StrikeSlipPotentials(s, z, beta, theta);

%%% Assemble.
Kernel =  Shear_1 - Shear_2;

    function Shear = StrikeSlipPotentials(s, z, beta, theta)
    
    %%% z_0 is the location of the dislocation. I.E. on-fault coordinate.
        z_0 = s*exp(1i*beta);
        zbar_0 = conj(z_0);
    
    %%% Potential.    
        Y = exp(1i*theta).*(1./(z - z_0) - 1./(z - zbar_0));
    
    %%% Shear stress due to shear dislocation
        Shear = imag(Y);
    end

end