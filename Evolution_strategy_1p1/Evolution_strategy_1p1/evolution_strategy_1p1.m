function [xk, k] = evolution_strategy_1p1(f, xkm1, sigma, ite, obj)
%%  [xk, k] = evolution_strategy_1p1(f, xkm1, sigma, ite)
%
%   This algorithm implements the simplest of all "evolution strategies"
%   (ES), the (1+1)-ES. In each iteration, one parent is used to create one
%   offspring using a Gaussian mutation operator (A random gaussian variable
%   with mean zero and standard deviation 'sigma').
%
%   INPUT DATA:
%
%   - f:     Function to minimize (handle function)
%   - xkm1:  Initial point (nx1 vector)
%   - sigma: Mutation strength (nx1 vector)
%   - ite:   Maximum number of iterations (positive integer number)
%   - obj:   Objective value
%
%   OUTPUT DATA:
%
%   - xk: Point which minimizes the function 'f' (nx1 vector)
%   - k:  Number of iterations (positive scalar)
%
%   Bibliography:
%
%   - KALYANMOY, Deb. "Multi-Objective optimization using evolutionary
%     algorithms". John Wiley & Sons, Ltd. Kanpur, India. 2001.
%
% -------------------------------------------------------
% | Developed by:   Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
%
%   Date: 20 - Sep - 2011

%% Beginning
n       = length(xkm1);           % 'n' states
xk      = zeros(n,ite);           % Pre-allocate space for 'xk'
maxE    = 1e-3;                   % maximum error
k       = 1;                      % counter
e       = 10;                     % initial high error
xk(:,1) = xkm1;                   % initial guessing
fx      = f(xk(:,1));             % evaluate function with 'xk'
while ((k < ite) && (e > maxE))
  y  = xk(:,k) + sigma.*randn(n,1); % create mutated solution 'y'
  fy = f(y);                        % evaluate function with mutated solution 'y'
  if (fy < fx)
    xk(:,k+1) = y;                  % update value of xkm1
    fx        = fy;                 % update value of f(xkm1)
  else
    xk(:,k+1) = xk(:,k);            % retain value of xkm1
  end
  e = abs(obj - fx);              % error
  k = k+1;                        % update counter
end

end
%% END