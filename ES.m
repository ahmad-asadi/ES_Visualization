function [F, error, iteration, X]=ES(mu , lambda, n_x , limits, iterationCount , errorThreshold , fitness, sigma, dims, selectionMode,successRateStrategy, successRate, strategy,handles)

    iteration = 0 ;
    error = 1 ;
    
    prev_error = 1 ;
    
    XMesh = limits(dims(1),1) :0.5: limits(dims(1),2) ;
    YMesh = limits(dims(2),1) :0.5: limits(dims(2),2) ;
    TMesh = limits(dims(2),1) :0.5: limits(dims(2),2) ;
    RMesh = limits(dims(2),1) :0.5: limits(dims(2),2) ;
    FMesh = limits(dims(2),1) :0.5: limits(dims(2),2) ;
    
    ZMesh = zeros(size(XMesh , 2) , size(YMesh , 2)) ;
    for i = 1 : size(XMesh , 2)
        for j = 1 : size(YMesh , 2)
            ZMesh(j,i) = fitness([XMesh(1,i) , YMesh(1,j)]) ;
        end
    end

    X = generate_population(mu, n_x , limits);
    F = evaluate_population_fitness(X , fitness) ;
    axes(handles.population) ;
        axis([limits(1,1) , limits(1,2) , limits(2,1) , limits(2,2) , 0 , 1])  ;
    rotate3d on ;
    while(error > errorThreshold && iteration < iterationCount)
        set(handles.currentIteration, 'string' , iteration) ;
        selected_parents = select_parents(X,F,mu,lambda) ;
        children = recombination(selected_parents, sigma , lambda) ;
        children = mutate(children);
        children_F = evaluate_population_fitness(children , fitness) ;
        [X,F] = select_population(selected_parents , F , children, children_F, mu , selectionMode) ;
        prev_error = [prev_error error] ;
        error = calculate_error(F) ;
        set(handles.currentError , 'string' , error) ;
        axes(handles.error) ;        
        if(iteration > 2)
            plot([iteration - 1 , iteration] , prev_error([iteration - 1 , iteration])) ;
            hold on ;
            drawnow ;
        end
        iteration = iteration + 1 ;
        axes(handles.population) ;
        surf(XMesh, YMesh, ZMesh) ;
        hold on;
        scatter3(X(:,dims(1)) , X(:,dims(2)) , F , [], repmat([1,0,0],size(X,1),1) , 'filled') ;
        drawnow;
        hold off; 
        axes(handles.contours) ;
        contour(XMesh , YMesh , ZMesh , 30);
        hold on ;
        scatter3(X(:,dims(1)) , X(:,dims(2)) , F , [], repmat([1,0,0],size(X,1),1) , 'filled') ;
        drawnow;
        hold off;
    end
    axes(handles.error) ;
    hold off ;
%     rotate3d on ;
%     if(iteration < iterationCount) % reached error threshold
%         iteration = iteration - 1 ;
%     end
%     plot(prev_error(max(1,end - iterationCount + 3) : end)) ;
    
    
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

function children = recombination(selected_parents, sigma, lambda)
    children = zeros(lambda , size(selected_parents,2)) ;
    for i = 1 : lambda
        children(i,:) = selected_parents(mod(i,size(selected_parents,1)) + 1 , : ) + sigma * rand(1 , size(selected_parents,2)) ;
    end
end

function children = mutate(children)
end

function [X,F] = select_population(parents, parents_fitness , children , children_fitness , mu , selectionMode) % selectionMode = 1 -> mu + lambda, selectionMode ~= 1 -> mu, lambda
    X = zeros(mu, size(parents,2)) ;
    F = zeros(mu, 1) ;
    if(selectionMode == 1)
        current_population = [parents;children] ;
        current_fitnesses = [parents_fitness; children_fitness] ;
    else
        current_population = children ;
        current_fitnesses = children_fitness ;
    end
    [~,sorted_fitnesses] = sort(current_fitnesses , 'descend') ;
    for i = 1 : mu
        X(i,:) = current_population(sorted_fitnesses(i),:) ;
        F(i,1) = current_fitnesses(sorted_fitnesses(i),1) ;
    end
end

function error = calculate_error(F)
    error = sum(F.^2) ;
end