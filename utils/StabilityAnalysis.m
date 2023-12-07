function Out_S = StabilityAnalysis(p, calc, geom)
%%% This is the dimensionless version.

%%% Geometry case.
if strcmp(geom, 'Thrust Fault')
    geometry = 'halfspace';
    sign = 1;
elseif strcmp(geom ,'Normal Fault')
    geometry = 'halfspace';
    sign = -1;
elseif strcmp(geom ,'Full Space')
    geometry = 'fullspace';
else
    geometry = 'layer';
end

%%% Geometrical parameters.
N = p.Geometry.GridPoints;
Xi = p.Geometry.FaultCoordinates;
dxi = p.Geometry.GridSpacing;

%%% Friction parameters.
a = p.Friction.a;
b = p.Friction.b;
mu = p.Friction.mu_ref;

%%% Length scales.
L_b = p.Friction.L_b;

switch geometry
%%% geometry-specific parameters. need different cases for thrust and
%%% normal faults.
%%% See notes on p.65 in Bienfang1 for half-space equations.    
    case 'halfspace'
        K_s = (1/(2*pi))*p.Geometry.ShearKernel;
        K_n = sign*(1/(2*pi))*p.Geometry.NormalKernel;

    case 'fullspace'
        K_s = (1/(2*pi))*FullSpace((1:N), (1:N), Xi');                       
        K_n = 0;

    case 'layer'
        BurialDepth = p.Geometry.BurialDepth;
        D2 = full(three_point_centered_uni_D2(Xi(1), Xi(N), N));
        D2(1,2) = 1/dxi^2;
        D2(N, N-1) = 1/dxi^2;
        K_s = 2*BurialDepth*D2;
        K_n = 0;
end

%%% Build stability eigenvalue problem matrix.
%%% Identity matrix.
I = eye(N);
%%% Upper left block.
A_ul = (b/a)*(L_b*(K_s - mu*K_n) + I);
%%% Upper right block.
A_ur = (b/a)*I;
%%% Lower left block.
A_ll = -I;
%%% Lower right block.
A_lr = -I;
%%% Construct marix.
A = [A_ul A_ur; A_ll A_lr];

switch calc
    case 'full'
%%% Get all the eigenvalues and eigenvectors.
        [V, D] = eig(A);
        lambda = diag(D);
%%% Compute determinant and trace of the Jacobian.
        Tr = trace(A);
        De = det(A);
%%% The eigenvalues are sets of complex conjugates. Look at the
%%% eigenvalues with positive imaginary parts, and exclude the eigenvalue
%%% lambda = 0. WHAT INFO IS BEING DISCARDED HERE?
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
    H = diff(Xi(kk));

%%% Get the along fault location and depths of the wavelengths. The
%%% location of each wavelength measurement is taken at the midpoint
%%% between successive peaks.
    Xi_kk = Xi(kk(1:end-1)) + H/2;

    switch geometry
        case 'halfspace'
%%% The depths are stored in p. Corresponding indices for Xi_kk would need
%%% to be computed. For now, just recompute the depths using the values of
%%% Xi_kk.
            beta = p.Geometry.FaultDip;
            alpha = p.Geometry.SurfaceSlope;
            Depth_kk = Xi_kk*(sin(pi*beta/180) + cos(pi*beta/180)*tan(pi*alpha/180));

        case 'fullspace'
            Depth_kk = nan*ones(numel(Xi_kk), 1);

        case 'layer'
            Depth_kk = BurialDepth*ones(numel(Xi_kk), 1);
            
        case 'strikeslip'
            beta = p.Geometry.FaultDip;
            alpha = p.Geometry.SurfaceSlope;
            Depth_kk = Xi_kk*(sin(pi*beta/180) + cos(pi*beta/180)*tan(pi*alpha/180));
    end
    
%%% Return the results.
    Out_S = struct('Case', geometry, 'EigenValues', lambda, 'EigenVectors', V,...
        'lambda_1', lambda_1, 'vector_1', vector_1, 'WaveLengths', H,...
        'WaveLengthLoc', [Xi_kk, Depth_kk], 'Stability', 0, 'Determinant',...
        De, 'Trace', Tr);
else
    Out_S = struct('Case', geometry, 'EigenValues', l_pos, 'EigenVectors', V_pos,...
        'Stability', 1, 'Determinant', De, 'Trace', Tr);
end

    case 'boundary'
%%% Find the eigenvalue with the largest real part.
        switch geometry            
            case 'fullspace'
                lambda = eigs(A, 1, 'largestreal');

            case 'halfspace'
                if N < 50
                    lambda = eigs(A, 1, 'largestreal');
                else
                    lambda = eigs(A, 1, 'largestreal', 'SubspaceDimension', 50,...
                        'Tolerance', 1e-11);
                end

            case 'layer'
                lambda = eigs(A, 1, 'largestreal', 'SubspaceDimension', 50,...
                    'Tolerance', 1e-11);
        end

        if real(lambda) > 0
            Out_S = struct('Case', geometry, 'EigenValue', lambda, 'Stability', 0);
        else
            Out_S = struct('Case', geometry, 'EigenValue', lambda, 'Stability', 1);
        end
end

%%% Full space shear stress operator.
function K = FullSpace(i, j, xi)
    dxi = mean(diff(xi));
    K = 1./(xi(j) - dxi/2 - xi(i)') - 1./(xi(j) + dxi/2 - xi(i)');
end

%%% Assign a phase angle to locations on the fault. 
% Xi_phase = nan(N,1);
% 
% for i = 1:numel(kk)-1
%     phase = linspace(0, 2*pi, numel(Xi(kk(i):kk(i+1))))';    
%     Xi_phase(kk(i):kk(i+1)-1) = phase(1:end-1);
% end
% 
% W = nan(N,1);
% for q = 1:numel(kk)-1
%     for i = kk(q):kk(q+1)-1
%         Xi_section1 = Xi(kk(q):kk(q+1)-1);
%         Xi_section2 = Xi(kk(q+1):kk(q+2)-1);
%         [~, j] = min(abs(Xi_phase(i) - Xi_phase(kk(q+1):kk(q+2)-1)));
%         W(i) = Xi_section2(j) - Xi_section1(i - kk(q) + 1);
%     end
% end


end
