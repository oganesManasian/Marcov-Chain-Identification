function X = MP_chain_1_v2(N_chain, Time, pi_a, x0)
    % Loading transition probabilities from chain 1, since we are using
    % them for base chain
    ld = load('P_hat_chain_1.mat');
    prob_matrix = ld.P_hat;
    
    X = zeros(Time + 1, N_chain);
    X(1, :) = x0; % Set initial state
    
    % Generating chain realisations
    prev_state = x0;
    for chain_ind = 1:N_chain
        for time = 2:Time + 1
            % Compute new state
            prev_state = X(time - 1, chain_ind);
            states_pair = chain_1(1, 2, prev_state);
            new_state = states_pair(2);
            
            % Compute acceptance
            nominator = pi_a(new_state) * prob_matrix(new_state, prev_state);
            denominator = pi_a(prev_state) * prob_matrix(prev_state, new_state);
            acceptance_rate = min(1, nominator / denominator);
            
            % Do or not do move
            do_move = randsrc(1, 1, [1, 0; acceptance_rate, 1 - acceptance_rate]);
            if do_move
                X(time, chain_ind) = new_state;
            else
                X(time, chain_ind) = prev_state;
            end
        end
    end
     
end