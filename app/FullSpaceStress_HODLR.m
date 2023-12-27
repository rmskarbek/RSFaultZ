function Kernel = FullSpaceStress_HODLR(i, j, xi, flag)

%%% Grid spacing.
    dxi = mean(diff(xi));

%%% Edge of fault element in negative xi-direction
    s_n = xi(j) - dxi/2;

%%% Edge of fault element in positive xi-direction.
    s_p = xi(j) + dxi/2;

%%% Assemble.
    if flag == 1
%%% Shear stress.
        Kernel = 1./(s_n - xi(i)') - 1./(s_p - xi(i)');
    else
%%% Normal stress. This should be the bimaterial operator.
        Kernel = 0;
    end

end