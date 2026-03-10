function [y, J] = damped_oscillator(theta, Set, t)
% DAMPED_OSCILLATOR_PHYSICAL
% Forward model and Jacobian for an underdamped mass–spring–damper system.
%
% INPUTS:
%   theta = [A; gamma; omega; phi]
%       A     : amplitude
%       gamma : damping coefficient
%       omega : angular frequency
%       phi   : phase
%   t     : column vector of time points (Nx1)
%
% OUTPUTS:
%   y : observations (Nx1)
%   J : Jacobian dy/dtheta (Nx4)

   % Unpack parameters
    A     = theta(1);
    gamma = theta(2);
    omega = theta(3);
    phi   = theta(4);

    % Ensure column vector
    t = t(:);

    % Common terms
    exp_term  = exp(-gamma * t);
    phase     = omega * t + phi;
    cos_term  = cos(phase);
    sin_term  = sin(phase);

    % Forward model
    y = A * exp_term .* cos_term;

    % Jacobian matrix
    % Each column corresponds to derivative w.r.t. a parameter
    J = zeros(length(t), 4);

    % dy/dA
    J(:,1) = exp_term .* cos_term;

    % dy/dgamma
    J(:,2) = -A * t .* exp_term .* cos_term;

    % dy/domega
    J(:,3) = -A * exp_term .* t .* sin_term;

    % dy/dphi
    J(:,4) = -A * exp_term .* sin_term;

    % clip
    y = y(Set);
    J = J(Set, :)'; 

    if any(isnan(y))
        error('nan encountered');
    elseif any(~isreal(y))
        error('imag encountered');
    end
end
