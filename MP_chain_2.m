function X = MP_chain_2(N_chain, Time, pi_a, x0)
    % Loading transition probabilities from chain 1, since we are using
    % them for base chain
    ld = load('P_hat_chain_2.mat');
    trans_prob = ld.P_hat;
    
    % Since state space is quite small (only 5 possible states, let's 
    % compute the whole new transition matrix
    state_size = size(trans_prob, 1);
    nominator = repmat(pi_a, state_size, 1) .* trans_prob';
    denominator = nominator';
    acceptance = min(ones(state_size) - eye(state_size), nominator ./ denominator);
    new_trans_prob = trans_prob .* acceptance;
    
    % Recomputing probabilities of self loops
    for i = 1:state_size
        new_trans_prob(i, i) = 1 - sum(new_trans_prob(i, :));
    end

    X = zeros(Time, N_chain);
    X(1, :) = x0; % Set initial state
    
    for time = 2:Time
        for state = 1:state_size
            % Find how many movements have to be done from particular state
            prev_states = X(time - 1, :);
            nb_needed_moves = sum(prev_states == state);
            
            % Generate movements according to new transtion probabilities
            movements = randsrc(1, nb_needed_moves, ...
                [1, 2, 3, 4, 5; new_trans_prob(state, :)]);
            
            % Fill X with new movements
            cur_states = X(time, :);
            cur_states(prev_states == state) = movements';
            X(time, :) = cur_states;
        end
    end
    
end