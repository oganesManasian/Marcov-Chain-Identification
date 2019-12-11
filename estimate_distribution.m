function pi_t = estimate_distribution(X, Time, state_size)
    % Estimates distribution of states for each time moment from 1 to Time
    % Returns matrix of shape (state_size, Time)

    assert(size(X, 1) == Time)
  
    pi_t = zeros(state_size, Time);
    n_chain = size(X, 2);
    
    for time = 1:Time
        % Count number of occurrences of particular state among all chains
        for i = 1:state_size
            pi_t(i, time) = sum(X(time, :) == i);
        end
        % Normalize by number of chains
        pi_t(:, time) = pi_t(:, time) / n_chain;
    end
end