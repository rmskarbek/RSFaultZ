function Out_S = StabilityWedge4(p)
%%% This version does not recompute the kernels.

%%% Geometry case.
if strcmp(p.Options.Geometry,'Thrust Fault')
    geom = 'halfspace';
else
    geom = 'fullspace';
end

%%% Geometrical parameters.
N = p.Geometry.GridPoints;
Xi = p.Geometry.FaultCoordinates;

%%% Material parameters.
v_0 = p.Material.PlateRate;
ShearMod = p.Material.ShearModulus;
Poisson = p.Material.PoissonRatio;
sigma_0 = p.Material.BackgroundStress;

%%% Friction parameters.
a = p.Friction.a;
b = p.Friction.b;
d_c = p.Friction.d_c;
mu = p.Friction.mu_ref;
 
%%% Identity and zero matrices.
I = eye(N);
Z = zeros(N);

switch geom
%%% See notes on p.65 in Bienfang1 for half-space equations.    
  case 'halfspace'
%%% geometry-specific parameters.
    alpha = p.Geometry.SurfaceSlope;
    beta = p.Geometry.FaultDip;
    K_s = (ShearMod/(2*pi*(1 - Poisson)))*p.Geometry.ShearKernel;
    K_n = (ShearMod/(2*pi*(1 - Poisson)))*p.Geometry.NormalKernel;

%%% Build stability eigenvalue problem matrix.
%%% Upper left block.
    A_ul = (v_0/(a*sigma_0))*K_s + (b*v_0/(a*d_c))*I...
        - (mu*v_0/(a*sigma_0))*K_n;
%%% Upper right block.
    A_ur = (b*v_0^3/(a*d_c^2))*I;
%%% Lower left block.
    A_ll = -(1/v_0)*I;
%%% Lower right block.
    A_lr = -(v_0/d_c)*I;
%%% Construct marix.
    A = [A_ul A_ur; A_ll A_lr];

%%% See notes on p.64 in Bienfang1 for full space equations.
  case 'fullspace'
%%% Geometry-specific parameters.
    K = (ShearMod/(2*pi*(1 - Poisson)))*FullSpace((1:N), (1:N), Xi');
%%% Build stability eigenvalue problem matrix.
%%% Upper left block.
    A_ul = (v_0./(a*sigma_0)).*K + (b*v_0./(a*d_c)).*I;
%%% Upper right block.
    A_ur = (b*v_0^3./(a*d_c^2)).*I;
%%% Lower left block.
    A_ll = -(1/v_0)*I;
%%% Lower right block.
    A_lr = -(v_0/d_c)*I;
%%% Construct marix.
    A = [A_ul A_ur; A_ll A_lr];

    case 'strikeslip'
%%% Geometry-specific parameters.
    alpha = p.Geometry.SurfaceSlope;
    beta = p.Geometry.FaultDip;
    K = (ShearMod/(2*pi*(1 - Poisson))).*p.Geometry.ShearKernel;
%%% Build stability eigenvalue problem matrix.
%%% Upper left block.
    A_ul = (v_0./(a*sigma_0)).*K + (b*v_0./(a*d_c)).*I;
%%% Upper right block.
    A_ur = (b*v_0^3./(a*d_c^2)).*I;
%%% Lower left block.
    A_ll = -(1/v_0)*I;
%%% Lower right block.
    A_lr = -(v_0/d_c)*I;
%%% Construct marix.
    A = [A_ul A_ur; A_ll A_lr];
end


%%% Get the eigenvalues and eigenvectors.
[V, D] = eig(A);
lambda = diag(D);

%%% eigs is about 5 times faster, but have to be sure that you get the
%%% smallest, positive real eigenvector. Maybe you can write own solver
%%% that does this.
% tic; [Vs, Ds] = eigs(A, 10, 'largestreal'); toc

%%% Remove the H_star condition, and replace with some sort of break
%%% statement that depends on if any eigenvalues have positive real parts.
%%% Output a message like "Fault is stable: no eigenvalues with positive
%%% real part".

%%% Find the nucleation lengths for VW fault.
% if p.Friction.H_star > 0
%%% The eigenvalues are sets of complex conjugates. Look at the
%%% eigenvalues with positive imaginary parts, and exclude the eigenvalue
%%% lambda = 0.

rr = imag(lambda) > 0 & abs(lambda) > 1e-15;
l_pos = lambda(rr);
%%% Get corresponding eigenvectors.
V_pos = V(:, rr);

%%% Determine if any eigenvalues have positive real parts.
if sum(real(l_pos) > 0) > 0

%%% Find the first eigenvalue and vector to cross the imaginary axis.
    jj = real(l_pos) == min(real(l_pos(real(l_pos)>0)));
    lambda_1 = l_pos(jj);
    vector_1 = V_pos(1:N, jj);

%%% The wavelength(s) of this eigenvector are the critical wavelength(s).
    [~, kk] = findpeaks(real(vector_1));
%    H = mean(diff(X(kk))); % this is from the fullspace stability code.
    H = diff(Xi(kk));

%%% Get the along fault location and depths of the wavelengths.
    Xi_kk = Xi(kk(1:end-1)) + H; % why is H added here? should it be H/2 to get the mid point between each peak in vector_1?
    switch geom
        case 'halfspace'
            Depth_kk = Xi_kk*(sin(pi*beta/180) + cos(pi*beta/180)*tan(pi*alpha/180));
        case 'fullspace'
            Depth_kk = nan*ones(numel(Xi_kk),1);
        case 'strikeslip'
            Depth_kk = Xi_kk*(sin(pi*beta/180) + cos(pi*beta/180)*tan(pi*alpha/180));
    end
    
%%% Return the results.
    Out_S = struct('Case', geom, 'EigenValues', lambda, 'EigenVectors', V,...
        'lambda_1', lambda_1, 'vector_1', vector_1, 'WaveLengths', H,...
        'WaveLengthLoc', [Xi_kk, Depth_kk], 'Stability', 0);
else
    Out_S = struct('Case', geom, 'EigenValues', lambda, 'EigenVectors', V,...
        'Stability', 1);
end

function K = FullSpace(i, j, xi)
    dxi = mean(diff(xi));
    K = 1./(xi(j) - dxi/2 - xi(i)') - 1./(xi(j) + dxi/2 - xi(i)');
end

end
