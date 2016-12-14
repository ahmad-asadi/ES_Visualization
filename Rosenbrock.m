function y = Rosenbrock(X)
    y = 0 ;
    for i = 1 : size(X,2) - 1
        y = y + 100 * (X(1,i+1) - X(1,i)^2)^2 + (1 - X(1,i))^2 ;
    end
    y = -y/100000 ;
end