%% EXAMPLE
clear, clc, close all

%%  Rosenbrock's function
f   = @(x,u) (1-x(1))^2 + 100*(x(2)-x(1)^2)^2;    % Minimum f(1,1) = 0
x0  = [3; -3];                                    % initial state
obj = 0;                                          % objective value

%% Allocate space in memory for variables
ite   = 1000;               % number of iterations
n_x   = length(x0);         % 'n_x' states
sigma = 0.5*ones(n_x,1);    % number of covariances

%% Compute minimum through 1p1-ES
[xk, k]  = evolution_strategy_1p1(f, x0, sigma, ite, obj);

%%  Plot results
x1 = zeros(1,k);
x2 = zeros(1,k);
for i = 1:k
  x1(i) = xk(1,i);
  x2(i) = xk(2,i);
end

figure
plot(x1)
xlabel('Iterations', 'FontSize', 16)
ylabel('x_{1}', 'FontSize', 16)
title('Convergence x_{1}', 'FontSize', 18)
grid on
figure
plot(x2)
xlabel('Iterations', 'FontSize', 16)
ylabel('x_{2}', 'FontSize', 16)
title('Convergence x_{2}', 'FontSize', 18)
grid on

%% END