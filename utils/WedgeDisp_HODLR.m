function Kernel = WedgeDisp_HODLR(i, j, xi, z, beta, alpha, G, nu, flag)
%%% The function computes the complex potentials for displacement
%%% calculations.

%%% Fault dip in z-plane.
    beta = pi - beta*pi/180;

%%% Wedge transform exponent.
    q = pi/(pi + alpha*pi/180);

%%% Dislocation strength, located on the fault at angle beta.
    gamma = -1i*G*exp(1i*beta)/(4*pi*(1 - nu));

%%% Build discretized elastic displacement operator assuming piecewise 
%%% continuous slip distribution.

%%% Grid spacing.
    dxi = mean(diff(xi));

%%% Up-dip edge of fault element.
    s = xi(i) - dxi/2;
    [u_x1, u_y1] = WedgePotentials(s, z(j), q, beta, gamma, nu, G);

%%% Down-dip edge of fault element.
    s = xi(i) + dxi/2;
    [u_x2, u_y2] = WedgePotentials(s, z(j), q, beta, gamma, nu, G);

%%% Assemble.
    if flag == 1
        Kernel = u_x1 - u_x2;
    else
        Kernel = u_y1 - u_y2;
    end

    function [u_x, u_y] = WedgePotentials(s, z, q, beta, gamma, nu, G)
        gammabar = conj(gamma);

    %%% z_0 is the location of the dislocation. I.E. on-fault coordinate.
        z_0 = s*exp(1i*beta);
        zbar_0 = conj(z_0);

    %%% Transformed half-space coordinates.
        zeta = z.^q;
        zeta_0 = z_0.^q;
        zetabar_0 = zbar_0.^q;

    %%% Omega (upper case).
        O = gamma*(bclog(zeta - zeta_0,beta) - bclog(zeta - zetabar_0,-beta))...
            + gammabar.*(zeta_0 - zeta)./(zeta - zetabar_0);
    %%% Omega prime.
        Op = gamma*(1./(zeta - zeta_0) - 1./(zeta - zetabar_0))...
            - gammabar*(1./(zeta - zetabar_0) + (zeta_0 - zeta)./(zeta - zetabar_0).^2);
    %%% omega (lower case).
        w = gammabar*(bclog(zeta - zeta_0,beta) - bclog(zeta - zetabar_0,-beta))...
            - gamma*zetabar_0./(zeta - zeta_0) + (gamma + gammabar)*zeta./(zeta - zetabar_0)...
            + gammabar*(zeta_0 - zeta).*zeta./(zeta - zbar_0).^2;

    %%% Complex displacement.
        u_c = (1./(2*G)).*((3 - 4*nu).*O - q*conj(zeta).*(z./conj(z)).*conj(Op)...
            - conj(w));
    %%% Displacement in positive x-direction.
        u_x = real(u_c);
    %%% Displacement in positive y-direction.
        u_y = imag(u_c);

    end

    function [u_x, u_y] = WedgePotentials2(s, z, q, beta, gamma, nu, G)
        gammabar = conj(gamma);

    %%% z_0 is the location of the dislocation. I.E. on-fault coordinate.
        z_0 = s*exp(1i*beta);
        zbar_0 = conj(z_0);

    %%% Transformed half-space coordinates.
        zeta = z.^q;
        zeta_0 = z_0.^q;
        zetabar_0 = zbar_0.^q;

    %%% Omega (upper case).
        O = gamma*(bclog(zeta - zeta_0,beta) - bclog(zeta - zetabar_0,-beta))...
            + gammabar.*(zeta_0 - q*zeta)./(zeta - zetabar_0);
    %%% Omega prime.
        Op = q*zeta.^(1-1/q).*(gamma*(1./(zeta - zeta_0) - 1./(zeta - zetabar_0))...
            + gammabar*(q*zetabar_0 - zeta_0)./(zeta - zetabar_0).^2);
    %%% omega (lower case).
        w = gammabar*(bclog(zeta - zeta_0,beta) - bclog(zeta - zetabar_0,-beta))...
            - gamma*zetabar_0./(zeta - zeta_0) + q*(gamma + q*gammabar)*zeta./(zeta - zetabar_0)...
            + q*gammabar*(zeta_0 - q*zeta).*zeta./(zeta - zbar_0).^2;

    %%% Complex displacement.
        u_c = (1./(2*G)).*((3 - 4*nu).*O - q*conj(zeta).*(z./conj(z)).*conj(Op)...
            - conj(w));
    %%% Displacement in positive x-direction.
        u_x = real(u_c);
    %%% Displacement in positive y-direction.
        u_y = imag(u_c);

    end

    function x = bclog(z, theta)
    %%% Zp locates the branch cut. Zp2 shifts the displacements.
    
    
    % arg[z_, σ_: - Pi] := Arg[z Exp[-I (σ + Pi)]] + σ + Pi;
    % log[z_, σ_: - Pi] := Log[Abs[z]] + I arg[z, σ]
    
        x = log(abs(z)) + 1i*angle(z.*exp(-1i*(theta + pi))) + 1i*pi;
    
    %    x = log(abs(z)) + 1i*angle(z.*exp(-1i*(theta + pi))) + theta + pi;% + 1i*pi;
    %    x = log(abs(z)) + 1i*angle(z.*exp(1i*(theta+pi))) - 1i*theta + 1i*pi;
    %    x = log(abs(z)) + 1i*(angle(z) + pi);
    end

end