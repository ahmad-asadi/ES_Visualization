function [limits, fitness] = getFitnessFunction(func_num , n_x)
    switch func_num
        case 1
            fitness = @Rosenbrock ;
            limits = repmat([-5 5] , n_x , 1) ;
    end
end