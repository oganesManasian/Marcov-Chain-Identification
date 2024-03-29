function pi_t = estimate_distribution_old(X, time,  state_size)
    % Please use estimate_distribution function
    pi_t = zeros(state_size, 1);
    n_chain = size(X, 2);

    for i = 1:state_size
        pi_t(i) = sum(X(time, :) == i);
    end
    pi_t = pi_t / n_chain;
end