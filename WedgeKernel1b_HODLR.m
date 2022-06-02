function Kernel = WedgeKernel1b_HODLR(i, j, xi, beta, alpha, flag)

I = sqrt(-1);
%%% Along-fault spatial coordinate.
beta = pi - beta*pi/180;
q = pi/(pi + alpha*pi/180);
%theta = 0;
%theta = pi + alpha*pi/180;
theta = beta;

gamma = -I*exp(I*beta)/2;
gammabar = I*exp(-I*beta)/2;

%%% Variables s are the integrands
dxi = mean(diff(xi));

%%% Build discretized elastic stress operator assuming piecewise continuous
%%% slip distribution.             
s = xi(j) - dxi/2;

z = xi(i)'*exp(I*theta);
%   zbar = Xi(i)*exp(-I*theta);
z_0 = s*exp(I*beta);
zbar_0 = s*exp(-I*beta);

%%% Potentials.
F = -q*z.^(q-1).*(gamma./(z.^q - z_0.^q)...
    - gamma./(z.^q - zbar_0.^q)...
    - q*gammabar.*(z_0 - zbar_0).*zbar_0.^(q-1)./(z.^q - zbar_0.^q).^2);

Fp = q*gamma.*z.^(q-2).*(z.^q + (q-1)*z_0.^q)./(z.^q - z_0.^q).^2 ...
    - q*gamma.*z.^(q-2).*(z.^q + (q-1)*zbar_0.^q)./(z.^q - zbar_0.^q).^2 ...
    - q^2*gammabar.*z.^(q-2).*zbar_0.^(q-1).*(z_0 - zbar_0)...
    .*((q+1)*z.^q + (q-1)*zbar_0.^q)./(z.^q - zbar_0.^q).^3;

Y = q*z.^(q-1).*((zbar_0.^q./(z.^q - zbar_0.^q).^2 - zbar_0.*z_0.^(q-1)./(z.^q - z_0.^q).^2).*q.*gamma...
    + q^2*gammabar.*zbar_0.^(q-1).*(zbar_0 - z_0).*(z.^q + zbar_0.^q)./(zbar_0.^q - z.^q).^3 ...
    + gammabar.*(1./(z.^q - zbar_0.^q) - 1./(z.^q - z_0.^q)));

%%% Shear stress due to shear dislocation
Shear_1 = -imag(F + conj(F) + exp(-2*I*theta).*(z.*conj(Fp) + conj(Y)));

%%% Normal stress due to shear dislocation. Sign is switched to account for
%%% solid mechanics sign conventions.
Normal_1 = -real(F + conj(F) + exp(-2*I*theta).*(z.*conj(Fp) + conj(Y)));

s = xi(j) + dxi/2;

z_0 = s*exp(I*beta);
zbar_0 = s*exp(-I*beta);

%%% Potentials.
F = -q*z.^(q-1).*(gamma./(z.^q - z_0.^q)...
    - gamma./(z.^q - zbar_0.^q)...
    - q*gammabar.*(z_0 - zbar_0).*zbar_0.^(q-1)./(z.^q - zbar_0.^q).^2);

Fp = q*gamma.*z.^(q-2).*(z.^q + (q-1)*z_0.^q)./(z.^q - z_0.^q).^2 ...
    - q*gamma.*z.^(q-2).*(z.^q + (q-1)*zbar_0.^q)./(z.^q - zbar_0.^q).^2 ...
    - q^2*gammabar.*z.^(q-2).*zbar_0.^(q-1).*(z_0 - zbar_0)...
    .*((q+1)*z.^q + (q-1)*zbar_0.^q)./(z.^q - zbar_0.^q).^3;

Y = q*z.^(q-1).*((zbar_0.^q./(z.^q - zbar_0.^q).^2 - zbar_0.*z_0.^(q-1)./(z.^q - z_0.^q).^2).*q.*gamma...
    + q^2*gammabar.*zbar_0.^(q-1).*(zbar_0 - z_0).*(z.^q + zbar_0.^q)./(zbar_0.^q - z.^q).^3 ...
    + gammabar.*(1./(z.^q - zbar_0.^q) - 1./(z.^q - z_0.^q)));

%%% Shear stress due to shear dislocation
Shear_2 = -imag(F + conj(F) + exp(-2*I*theta).*(z.*conj(Fp) + conj(Y)));

%%% Normal stress due to shear dislocation. Sign is switched to account for
%%% solid mechanics sign conventions.
Normal_2 = -real(F + conj(F) + exp(-2*I*theta).*(z.*conj(Fp) + conj(Y)));

%%% Assemble.
if flag == 1
    Kernel =  Shear_1 - Shear_2;
else
    Kernel = Normal_1 - Normal_2;
end