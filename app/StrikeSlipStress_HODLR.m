function Kernel = StrikeSlipStress_HODLR(i, j, xi, beta, alpha, flag)
%%% This is in my coordinates.
I = sqrt(-1);
%%% Along-fault spatial coordinate.
beta = (pi/180)*beta;
%theta = 0;
%theta = pi + alpha*pi/180;
theta = beta - pi/2;

%%% Variables s are the integrands
dxi = mean(diff(xi));

%%% Build discretized elastic stress operator assuming piecewise continuous
%%% slip distribution.             
s = xi(j) - dxi/2;

z = xi(i)'*exp(I*beta);
%   zbar = Xi(i)*exp(-I*theta);
z_0 = s*exp(I*beta);
zbar_0 = s*exp(-I*beta);

%%% Potential.
Y = exp(I*theta).*(1./(z - z_0) - 1./(z - zbar_0));

%%% Shear stress due to screw dislocation
Shear_1 = imag(Y);


%%% Edge 2.
s = xi(j) + dxi/2;

z_0 = s*exp(I*beta);
zbar_0 = s*exp(-I*beta);

%%% Potentials.
Y = exp(I*theta).*(1./(z - z_0) - 1./(z - zbar_0));

%%% Shear stress due to screw dislocation
Shear_2 = imag(Y);


%%% Assemble.
if flag == 1
    Kernel =  Shear_1 - Shear_2;
else
    Kernel = Normal_1 - Normal_2;
end
