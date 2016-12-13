function f = drawBenchmarkContours(func_num , X1, X2)
    [f, n_x , limits , obj] = get_benchmark_function(func_num) ;
    fprintf('Drawing function contours\n') ;
    fprintf('Number of dimensions: %d \n' , n_x);
    fprintf('Objective value: %d\n' , obj) ;
    fprintf('Boundaries: ') ;
    disp(limits);
    
    if(n_x <= 2)
        X1 = 1 ;
        X2 = 2 ;
    end
    
    X1_limits = limits(X1,:) ;
    if(n_x > 1)
        X2_limits = limits(X2,:) ;
    end
    
    if(n_x > 1)
        ezsurf(f , X1_limits, X2_limits);
    else
        ezplot(f, X1_limits);
    end





function [f , n_x, limits, obj] = get_benchmark_function(fun)
switch fun
  case 1
    %%  Rosenbrock's function
    %   Minimum: f(1,1) = 0
    f      = @(x,u) (1-x).^2 + 100*(u-x.^2).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 2
    %%  Himmelblau's function
    %   Minimum:  * f( 3,         2       ) = 0
    %             * f(-2.805118,  3.131312) = 0
    %             * f(-3.779310, -3.283186) = 0
    %             * f( 3.584428, -1.848126) = 0
    f      = @(x,u) (x.^2 + u - 11).^2 + (x + u.^2 - 7).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
       
  case 3
    %%  Rastrigin function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 20 + (x.^2 + u.^2) - 10*(cos(2*pi*x) + cos(2*pi*u));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 4
    %%  Sphere function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) x.^2 + u.^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
    
  case 5
    %%  Ackley's function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) -20*exp(-0.2*sqrt(0.5*(x.^2 + u.^2))) - exp(0.5*(cos(2*pi*x) + cos(2*pi*u))) + 20 + exp(1);
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
  
  case 6
    %%  F1 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0) = 1
    f      = @(x,u) abs(x) + cos(x);
    n_x    = 1;                           % 'n_x' states
    limits = repmat([-20 20], n_x, 1);    % Boundaries
    obj    = 1;                           % objective value (f(x_min) = obj)
    
  case 7
    %%  F2 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0) = 0
    f      = @(x,u) abs(x) + sin(x);
    n_x    = 1;                           % 'n_x' states
    limits = repmat([-20 20], n_x, 1);    % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 8
    %%  F5 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(zeros(1,n_x)) = -20
    n_x    = 2;                           % 'n_x' states
    f      = @(x,u) sum(abs(x(1:n_x,:)) - 10*cos(sqrt(abs(10*x(1:n_x,:)))));
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = -20;                         % objective value (f(x_min) = obj)

  case 9
    %%  F6 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(9.6204) = -100.22
    f      = @(x,u) (x.^2 + x).*cos(x);
    n_x    = 1;                           % 'n_x' states
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = -100.22;                     % objective value (f(x_min) = obj)
    
  case 10
    %%  F7 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(9.039,8.668) = -18.5547
    f = @(x,u) x.*sin(4*x) + 1.1*u.*sin(2*u);
    n_x    = 2;                           % 'n_x' states
    limits = repmat([0 10], n_x, 1);      % Boundaries
    obj    = -18.5547;                    % objective value (f(x_min) = obj)

  case 11
    %%  F8 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(9.039,8.668) = -18.5547
    f = @(x,u) u.*sin(4*x) + 1.1*x.*sin(2*u);
    n_x    = 2;                           % 'n_x' states
    limits = repmat([0 10], n_x, 1);      % Boundaries
    obj    = -18.5547;                    % objective value (f(x_min) = obj)
  
  case 12
    %%  F11 (See Appendix I (Haupt, 2004)).
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 1 + ((x.^2 + u.^2)/4000) - cos(x).*cos(u);
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 13
    %%  F12 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(zeros(1,n_x)) = 0
    n_x    = 2;                           % 'n_x' states
    f      = @(x,u) 0.5 + ((sin(sqrt(sum(x(1:n_x,:).^2)))).^2 - 0.5)./(1 + 0.1*(sum(x(1:n_x,:).^2)));
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 14
    %%  F13 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0,0) = 0
    f      = @(x,u) (x.^2 + u.^2).^(0.25).*sin(30*((x + 0.5).^2 + u.^2).^0.1) + abs(x) + abs(u);
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-10 10], n_x, 1);    % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 15
    %%  F15 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(0,0.7688) = -345.36
    f      = @(x,u) -exp(-0.2*sqrt(x.^2 + u.^2) + 3*(cos(2*x) + sin(2*u)));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = -345.36;                     % objective value (f(x_min) = obj)

  case 16
    %%  F16 (See Appendix I (Haupt, 2004)).
    %   Minimum: f(-20,3.15) = -23.59
    f      = @(x,u) -x.*sin(sqrt(abs(x - (u+9)))) - (u+9).*sin(sqrt(abs(u + 0.5*x + 9)));
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-20 20], n_x, 1);    % Boundaries
    obj    = -23.59;                      % objective value (f(x_min) = obj)

  case 17
    %%  Other function (no reference)
    %   Minimum: f(0,0) = -20.09
    f      = @(x,u) -exp(-0.2*sqrt(x.^2 + u.^2) + 3*(cosd(2*x) + sind(2*u)));
    n_x    = 2;                         % 'n_x' states
    limits = repmat([-5 5], n_x, 1);    % Boundaries
    obj    = -20.09;                    % objective value (f(x_min) = obj)

  case 18
    %%  Other function (No reference)
    %   Minimum: f(zeros(1,n_x)) = 0
    n_x    = 2;                           % 'n_x' states
    f      = @(x,u) 0.5 + ((sqrt(sum(x(1:n_x,:).^2))).^2 - 0.5)./(1 + 0.1*(sum(x(1:n_x,:).^2)));
    limits = repmat([-5 5], n_x, 1);      % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)

  case 19
    %%  Beale's function
    %   Minimum: f(3,0.5) = 0
    f      = @(x,u) (1.5 - x + x.*u).^2 + (2.25 - x + x.*u.^2).^2 + (2.625 - x + x.*u.^3).^2;
    n_x    = 2;                           % 'n_x' states
    limits = repmat([-4.5 4.5], n_x, 1);  % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
   
  case 20
    %%  Rosenbrock's function (Generalization to n_x dimensions)
    %   Minimum: f(1,1,...,1) = 0
    n_x    = 4;                           % 'n_x' states (variable)
    f      = @(x,u) sum(100*(x(2:2:n_x, :) - x(1:2:n_x-1, :).^2).^2 + (1 - x(1:2:n_x-1, :)).^2, 1);
    limits = repmat([-100 100], n_x, 1);  % Boundaries
    obj    = 0;                           % objective value (f(x_min) = obj)
    
  case 21
    %%  Bird's function
    %   Minimum:  * f( 4.701055751981055,  3.152946019601391) = -106.764537
    %             * f(-1.582142172055011, -3.130246799635430) = -106.764537
    f      = @(x,u) sin(x).*exp((1-cos(u)).^2) + cos(u).*exp((1-sin(x)).^2) + (x-u).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-2*pi 2*pi], n_x, 1);  % Boundaries
    obj    = -106.764537;                   % objective value (f(x_min) = obj)

  case 22
    %%  Booth's function
    %   Minimum: f(1,3) = 0
    f      = @(x,u) (x + 2*u - 7).^2 + (2*x + u - 5).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 23
    %%  Bukin4 function
    %   Minimum: f(-10,0) = 0
    f      = @(x,u) 100*u.^2 + 0.01*abs(x + 10);
    n_x    = 2;                             % 'n_x' states
    limits = [...
                  -15 -5
                  -3   3
             ];                             % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 24
    %%  Bukin6 function
    %   Minimum: f(-10,1) = 0
    f      = @(x,u) 100*sqrt(abs(u - 0.01*x.^2)) + 0.01*abs(x + 10);
    n_x    = 2;                             % 'n_x' states
    limits = [...
                  -15 -5
                  -3   3
             ];                             % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 25
    %%  Carrom Table function
    %   Minimum:  * f( 9.646157266348881,  9.646134286497169) = 24.1568155
    %             * f(-9.646157266348881,  9.646134286497169) = 24.1568155
    %             * f( 9.646157266348881, -9.646134286497169) = 24.1568155
    %             * f(-9.646157266348881, -9.646134286497169) = 24.1568155
    f      = @(x,u) -((cos(x).*cos(u).*exp(abs(1 - sqrt(x.^2 + u.^2)/pi))).^2)/30;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 24.1568155;                    % objective value (f(x_min) = obj)

  case 26
    %%  Chichinadze function
    %   Minimum: f(5.90133, 0.5) = -43.3159
    f      = @(x,u) x.^2 - 12*x + 11 + 10*cos(pi*x/2) + 8*sin(5*pi*x/2) - 1/sqrt(5)*exp(-((u - 0.5).^2)/2);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-30 30], n_x, 1);      % Boundaries
    obj    = -43.3159;                      % objective value (f(x_min) = obj)
    
  case 27
    %%  Cross function
    %   Minimum: f(NaN, NaN) = 0  ---> Too complicated to determine the real solution
    f      = @(x,u) (abs(sin(x).*sin(u).*exp(abs(100 - sqrt(x.^2 + u.^2)/pi))) + 1).^(-0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)
    
  case 28
    %%  Cross-in-tray function
    %   Minimum:  * f( 1.34940668535334, -1.34940668535334) = -2.062611870822739
    %             * f( 1.34940668535334,  1.34940668535334) = -2.062611870822739
    %             * f(-1.34940668535334,  1.34940668535334) = -2.062611870822739
    %             * f(-1.34940668535334, -1.34940668535334) = -2.062611870822739
    f      = @(x,u) -0.0001*(abs(sin(x).*sin(u).*exp(abs(100 - sqrt(x.^2 + u.^2)/pi))) + 1).^(0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -2.062611870822739;            % objective value (f(x_min) = obj)
    
  case 29
    %%  Cross-leg-table function
    %   Minimum:  f(NaN,NaN) = -1  ---> Too complicated to determine the real solution
    f      = @(x,u) -(abs(sin(x).*sin(u).*exp(abs(100 - sqrt(x.^2 + u.^2)/pi))) + 1).^(-0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -1;                            % objective value (f(x_min) = obj)

  case 30
    %%  Crowned cross function
    %   Minimum:  f(NaN,NaN) = 0.0001 ---> Too complicated to determine the real solution
    f      = @(x,u) 0.0001*(abs(sin(x).*sin(u).*exp(abs(100 - sqrt(x.^2 + u.^2)/pi))) + 1).^(0.1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -1;                            % objective value (f(x_min) = obj)
    
  case 31
    %%  Cube function
    %   Minimum:  f(1,1,...,1) = 0
    n_x    = 4;                             % 'n_x' states
    f      = @(x,u) sum(100*(x(2:n_x,:) - x(1:(n_x-1),:).^3).^2 + (1 - x(1:(n_x-1),:)).^2,1);
    limits = repmat([-20 20], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)
    
  case 32
    %%  Easom function
    %   Minimum:  f(pi,pi) = -1
    f      = @(x,u) -cos(x).*cos(u).*exp(-((x-pi).^2 + (u-pi).^2));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = -1;                            % objective value (f(x_min) = obj)
    
  case 33
    %%  Eggholder function (Just 2D)
    %   Minimum:  f(512, 404.2319) = -959.6407
    f      = @(x,u) sum(-(u+47).*sin(sqrt(abs(u+x/2+47)))-x.*sin(sqrt(abs(x-(u+47)))), 1);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-512 512], n_x, 1);    % Boundaries
    obj    = -959.6407;                     % objective value (f(x_min) = obj)
    
  case 34
    %%  Giunta function (Just 2D)
    %   Minimum:  f(0.45834282, 0.45834282) = 0.0602472184
    n_x    = 2;                             % 'n_x' states
    f      = @(x,u) 0.6 + sum(sin(16*x(:,:)/15 - 1) + sin(16*x(:,:)/15 - 1).^2 + sin(4*(16*x(:,:)/15 - 1))/50, 1);
    limits = repmat([-1 1], n_x, 1);        % Boundaries
    obj    = 0.0602472184;                  % objective value (f(x_min) = obj)

  case 35
    %%  Goldstein Price function
    %   Minimum:  f(0,-1) = 3
    f      = @(x,u) (1  + (x + u + 1).^2.*(19 - 14*x +  3*x.^2 - 14*u +  6*x.*u +  3*u.^2)).*...
                    (30 + (2*x - 3*u).^2.*(18 - 32*x + 12*x.^2 + 48*u - 36*x.*u + 27*u.^2));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-2 2], n_x, 1);        % Boundaries
    obj    = 3;                             % objective value (f(x_min) = obj)

  case 36
    %%  Griewank function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) (x.^2 + u.^2)/200 - cos(x).*cos(u/sqrt(2)) + 1;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 37
    %%  Helical Valley function
    %   Minimum:  f(1,0,0) = 0
    f      = @(x,u) 100*((x(3,:) - 10*atan2(u, x)/2/pi).^2 + (sqrt(x.^2 + u.^2) - 1).^2) + x(3,:).^2;
    n_x    = 3;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 38
    %%  Holder table function
    %   Minimum:  * f( 8.055023472141116, 9.664590028909654) = -19.2085026
    %             * f(-8.055023472141116, 9.664590028909654) = -19.2085026
    %             * f( 8.055023472141116,-9.664590028909654) = -19.2085026
    %             * f(-8.055023472141116,-9.664590028909654) = -19.2085026
    f      = @(x,u) -abs(sin(x).*cos(u).*exp(abs(1 - sqrt(x.^2 + u.^2)/pi)));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -19.2085026;                   % objective value (f(x_min) = obj)

  case 39
    %%  Levi13 function
    %   Minimum:  f(1,1) = 0
    f      = @(x,u) sin(3*pi*x).^2 + (x-1).^2.*(1 + sin(3*pi*u).^2) + (u-1).^2.*(1 + sin(2*pi*u).^2);
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 40
    %%  Matyas function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 0.26*(x.^2 + u.^2) - 0.48*x.*u;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 41
    %%  McCormick function
    %   Minimum:  f(-0.54719, -1.54719) = -1.9133
    f      = @(x,u) sin(x + u) + (x-u).^2 - 1.5*x + 2.5*u + 1;
    n_x    = 2;                             % 'n_x' states
    limits = [...
                -1.5   4
                -3     4
             ];                             % Boundaries
    obj    = -1.9133;                       % objective value (f(x_min) = obj)

  case 42
    %%  Schaffer1 function
    %   Minimum:  f(0, 0) = 0
    f      = @(x,u) 0.5  + (sin(x.^2 + u.^2).^2 - 0.5) ./ (1+0.001*(x.^2 + u.^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 43
    %%  Schaffer2 function
    %   Minimum:  f(0, 0) = 0
    f      = @(x,u) 0.5  + (sin(x.^2 - u.^2).^2 - 0.5) ./ (1+0.001*(x.^2 + u.^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 44
    %%  Schaffer3 function
    %   Minimum:  f(0, 1.25313) = 0.00156685
    f      = @(x,u) 0.5  + (sin(cos(abs(x.^2 - u.^2))).^2 - 0.5) ./ (1+0.001*(x.^2 + u.^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0.00156685;                    % objective value (f(x_min) = obj)

  case 45
    %%  Schaffer4 function
    %   Minimum:  f(0, 1.25313) = 0.292579
    f      = @(x,u) 0.5  + (cos(sin(abs(x.^2 - u.^2))).^2 - 0.5) ./ (1+0.001*(x.^2 + u.^2)).^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0.292579;                      % objective value (f(x_min) = obj)

  case 46
    %%  Pen-Holder function
    %   Minimum:  * f(-9.64617, 9.64617) = -0.96353
    %             * f( 9.64617,-9.64617) = -0.96353
    %             * f( 9.64617, 9.64617) = -0.96353
    %             * f(-9.64617,-9.64617) = -0.96353
    f      = @(x,u) -exp(-(abs(cos(x).*cos(u).*exp(abs(1 - sqrt(x.^2 + u.^2)/pi)))).^(-1));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-11 11], n_x, 1);      % Boundaries
    obj    = -0.96353;                      % objective value (f(x_min) = obj)

  case 47
    %%  Powell function
    %   Minimum:  f(0,0,0,0) = 0
    f      = @(x,u) (x + 10*u).^2 + 5*(x(3,:) - x(4,:)).^2 + (u - 2*x(3,:)).^4 + 10*(x - x(4,:)).^4;
    n_x    = 4;                             % 'n_x' states
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 48
    %%  Schweffel function
    %   Minimum:  f(420.9687, 420.9687) = -837.9658
    f      = @(x,u) -x.*sin(sqrt(abs(x))) -u.*sin(sqrt(abs(u)));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-500 500], n_x, 1);    % Boundaries
    obj    = -837.9658;                     % objective value (f(x_min) = obj)

  case 49
    %%  Sine Envelope Sine function
    %   Minimum:  f(0,0,...,0) = 0
    n_x    = 4;                             % 'n_x' states (variable)
    f      = @(x,u) sum((sin(sqrt(x(1:n_x-1, :).^2 + x(2:n_x, :).^2)).^2 - 0.5)./(1 + 0.001*(x(1:n_x-1, :).^2 + x(2:n_x, :).^2)).^2 + 0.5, 1);
    limits = repmat([-100 100], n_x, 1);    % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 50
    %%  Six Hump Camel Back function
    %   Minimum:  * f( 8.984201368301331e-2,-7.126564032704135e-001) = -1.031628453489877
    %             * f(-8.984201368301331e-2, 7.126564032704135e-001) = -1.031628453489877
    f      = @(x,u) (4 - 2.1*x.^2 + x.^4/3).*x.^2 + x.*u + (4*u.^2 - 4).*u.^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = -1.031628453489877;            % objective value (f(x_min) = obj)

  case 51
    %%  Styblinski-Tang function
    %   Minimum:  f(-2.903534, -2.903534, ..., -2.903534) = -39.16599*n_x
    n_x    = 2;                             % 'n_x' states (Variable)
    f      = @(x,u) sum(x(:,:).^4 - 16*x(:,:).^2 + 5*x(:,:), 1)/2;
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = -39.16599 * n_x;               % objective value (f(x_min) = obj)

  case 52
    %%  Test tube Holder function
    %   Minimum:  f(pi/2,0) = -10.8723
    f      = @(x,u) -4*abs(sin(x).*cos(u).*exp(abs(cos((x.^2 + u.^2)/200))));
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-10 10], n_x, 1);      % Boundaries
    obj    = -10.872299901558;              % objective value (f(x_min) = obj)

  case 53
    %%  Three Hump Camel function
    %   Minimum:  f(0,0) = 0
    f      = @(x,u) 2*x.^2 - 1.05*x.^4 + x.^6/6 + x.*u + u.^2;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 54
    %%  Wood function
    %   Minimum:  f(1,1,1,1) = 0
    f      = @(x,u) 100*(x.^2 - u).^2 + (x - 1).^2 + (x(3,:) - 1).^2 + 90*(x(3,:).^2 - x(4,:)).^2 + ...
                    10.1*((u - 1).^2 + (x(4,:) - 1).^2)  + 19.8*(u - 1).*(x(4,:) - 1);
    n_x    = 4;                             % 'n_x' states
    limits = repmat([-50 50], n_x, 1);      % Boundaries
    obj    = 0;                             % objective value (f(x_min) = obj)

  case 55
    %%  Zetti function
    %   Minimum:  f(-0.0299, 0) = -0.003791
    f      = @(x,u) (x.^2 + u.^2 - 2*x).^2 + x/4;
    n_x    = 2;                             % 'n_x' states
    limits = repmat([-5 5], n_x, 1);        % Boundaries
    obj    = -0.003791;                     % objective value (f(x_min) = obj)

  otherwise
    error('Not supported equation');
end
