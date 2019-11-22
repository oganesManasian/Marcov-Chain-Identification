function P = estimate_transition_matrix(X, time, state_size)
    % Computes transition rate for each possible path at particular time 
    P = zeros(state_size, state_size);
    n_chain = size(X, 2);
    
    for chain = 1:n_chain - 1
        i = X(time, chain);
        j = X(time + 1, chain);
        P(i, j) = P(i, j) + 1;
    end
    P = P / n_chain;
end