function P = estimate_transition_matrix(X, time, state_size)
    % Computes transition rate for each possible path at particular time 
    % Returns transition matrix if shape (state_size, state_size)
    
    P = zeros(state_size, state_size);
    n_chain = size(X, 2);
    
    % Count number of transitions for each possible path
    for chain = 1:n_chain - 1
        i = X(time, chain);
        j = X(time + 1, chain);
        P(i, j) = P(i, j) + 1;
    end
    
    % Normalise row wise
    for i = 1:state_size
        P(i, :) = P(i, :) / sum(P(i, :));
    end
end