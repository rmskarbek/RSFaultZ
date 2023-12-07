function Kernel = HalfSpaceSmallAngle_HODLR(i, j, xi, beta, flag)
%%% Fault dip in z-plane.
    beta = beta*pi/180;

%%% Variables xi(j) are the integrands
    dxi = mean(diff(xi));

%%% Build discretized elastic stress operator assuming piecewise continuous
%%% slip distribution.
    z = xi(i)';
%%% Up-dip edge of fault element.
    s = xi(j) - dxi/2;
    K_s1 = beta.^2*s.*(z + s)./(z - s).^3;
    K_n1 = beta.^3.*s.^2.*(2*z + s)./(z - s).^4;

%%% Down-dip edge of fault element.
    s = xi(j) + dxi/2;
    K_s2 = beta.^2*s.*(z + s)./(z - s).^3;
    K_n2 = beta.^3.*s.^2.*(2*z + s)./(z - s).^4;

%%% Assemble.
    if flag == 1
        Kernel = K_s1 - K_s2;
    else
        Kernel = K_n1 - K_n2;
    end