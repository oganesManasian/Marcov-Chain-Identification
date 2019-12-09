function X = MP_chain_1(N_chain, Time, pi_a, x0)
    % Loading transition probabilities from chain 1, since we are using
    % them for base chain
    ld = load('P_hat_chain_1.mat');
    prob_matrix = ld.P_hat';
    
    % Computing new transition matrix
    state_size = size(prob_matrix, 1);
    nominator = repmat(pi_a, state_size, 1) .* prob_matrix';
    denominator = nominator';
    acceptance = min(ones(state_size) - eye(state_size), nominator ./ denominator);
    prob_matrix_new = prob_matrix .* acceptance;
    
    % Recomputing probabilities of self loops
    for i = 1:state_size
        prob_matrix_new(i, i) = 1 - sum(prob_matrix_new(i, :));
    end
    
	prob_matrix_new;
    
    % Generating chain realisations
    X = zeros(Time + 1, N_chain);
    X(1, :) = x0; % Set initial state
    
    for chain_ind = 1:N_chain
        for time = 2:Time + 1
            prev_state = X(time - 1, chain_ind);
            trans_prob_from_prev_state = prob_matrix_new(prev_state, :);
            X(time, chain_ind) = randsrc(1, 1, [1, 2, 3, 4, 5; trans_prob_from_prev_state]);
        end
    end
     
end