function [F, error, iteration, X]=ES(mu , lambda, n_x , limits, iterationCount , errorThreshold , fitness, sigma)

    iteration = 0 ;
    error = Inf(1) ;
    
    X = generate_population(mu, n_x , limits);
    F = evaluate_population_fitness(X , fitness) ;
    
    while(error > errorThreshold && iteration < iterationCount)
        selected_parents = select_parents(X,F,mu,lambda) ;
        children = recombination(selected_parents, sigma) ;
        children = mutate(children);
        children_F = evaluate_population_fitness(children , fitness) ;
        [X,F] = select_population(selected_parents , F , children, children_F, lambda) ;
        prev_error = error ;
        error = calculate_error(F) ;
%         if(iteration > 2)
%             plot([iteration - 1 , iteration] , [prev_error, error]) ;
%             hold on ;
%             drawnow ;
%         end
        iteration = iteration + 1 ;
        ezsurf('1/(1+x1^2 + x2^2)' , [limits(1,1) , limits(1,2) , limits(2,1) , limits(2,2)]) ;
        hold on ;
        scatter(X(:,1) , X(:,2)) ;
        axis([limits(1,1) , limits(1,2) , limits(2,1) , limits(2,2)])  ;
        drawnow;
%        hold off;
    end
    
end

function X = generate_population(mu , n_x , limits)
    X = zeros(mu, n_x) ;
    for i = 1 : n_x
        X(:,i) = rand(mu,1) * (limits(i,2) - limits(i,1)) + limits(i,1) ;
    end
end

function F = evaluate_population_fitness(X , fitness) 
    F = zeros(size(X,1) , 1) ;
    for i = 1:size(X,1)
        F(i) = fitness(X(i,:)) ;
    end
end

function parents = select_parents(X , F , mu , lambda) 
    parents = X ;
end

function children = recombination(selected_parents, sigma)
    children = selected_parents + sigma * rand(size(selected_parents,1) , size(selected_parents,2)) ;
end

function children = mutate(children)
end

function [X,F] = select_population(parents, parents_fitness , children , children_fitness , lambda)
    X = zeros(lambda, size(parents,2)) ;
    F = zeros(lambda, 1) ;
    current_population = [parents;children] ;
    current_fitnesses = [parents_fitness; children_fitness] ;
    [~,sorted_fitnesses] = sort(current_fitnesses , 'descend') ;
    for i = 1 : lambda
        X(i,:) = current_population(sorted_fitnesses(i),:) ;
        F(i,1) = current_fitnesses(sorted_fitnesses(i),1) ;
    end
end

function error = calculate_error(F)
    error = abs(1 - F(1)) ;
end