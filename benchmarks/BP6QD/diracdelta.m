function d = diracdelta(x, h)
%%% Divyaprakash (2024). Smoothed Dirac Delta Function.
%%% (https://www.mathworks.com/matlabcentral/fileexchange/157136-smoothed-dirac-delta-function),
%%%% MATLAB Central File Exchange. Retrieved April 8, 2024.

% diracdelta - Compute the Dirac delta function at point x with mesh width h.
%   d = diracdelta(x, h) computes the Dirac delta function at the point x
%   with a specified mesh width h. The function is designed for use in
%   numerical simulations where h represents the grid spacing of a mesh. It
%   is defined such that it is differentiable and approaches an impulse as
%   h approaches zero.
%
% Input:
%   - x: Point(s) at which to evaluate the Dirac delta function.
%   - h: Mesh width (grid spacing) of the simulation grid.
%
% Output:
%   - d: Value of the Dirac delta function at point x with mesh width h.
% Author: Divyaprakash
%         Mechanical Engineer
% e-mail: divyaprakash.poddar@gmail.com
% Date  : 4 January 2024
    % Initialize dirac variable
    dirac = 1;
    % Iterate over dimensions (assuming x is a vector of length 2)
    for ii = 1:1
        % Compute normalized distance from x to the center of the distribution
        r = x(ii) / h;
        % Evaluate the piecewise function for the Dirac delta
        if (abs(r) <= 1.0)
            phi = (3.0 - 2 * abs(r) + sqrt(1.0 + 4.0 * abs(r) - 4.0 * r^2)) / 8.0;
        elseif and((abs(r) > 1.0), (abs(r) <= 2.0))
            phi = (5.0 - 2 * abs(r) - sqrt(-7.0 + 12.0 * abs(r) - 4.0 * r^2)) / 8.0;
        else
            phi = 0.0;
        end
        % Update dirac variable
        dirac = dirac * phi;
    end
    % Compute the final Dirac delta value by scaling and normalizing
    d = (1.0 / h^2) * dirac;
end
