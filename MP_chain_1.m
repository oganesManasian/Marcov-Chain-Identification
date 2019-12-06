function X = MP_chain_1(N_chain, Time, pi_a, x0)
    % Loading transition probabilities from chain 1, since we are using it
    % as base chain
    % TODO change to loading from file P_hat_chain_i.mat when 1C will be
    % done
    State_size = 5;
    X = chain_1(100000, 10, ones(1,5)/5);
    prob_matrix = estimate_transition_matrix(X, 1, State_size);
    
    % Computing new transition matrix
    matrix_size = size(prob_matrix, 1);
    nominator = repmat(pi_a, matrix_size, 1) .* prob_matrix';
    denominator = nominator';
    acceptance = min(ones(matrix_size) - eye(matrix_size), nominator ./ denominator);
    prob_matrix_new = prob_matrix .* acceptance;
    % Recomputing probabilities of self loops
    for i = 1:matrix_size
        prob_matrix_new(i, i) = 1 - sum(prob_matrix_new(i, :));
    end

    % Generating chain realisations
    X = zeros(N_chain, Time + 1);
    X(:, 1) = x0; % Set initial state
    for time = 2:Time + 1
        for chain_ind = 1:N_chain
            prev_state = X(chain_ind, time - 1);
            trans_prob_from_prev_state = prob_matrix_new(prev_state, :);
            X(chain_ind, time) = randsrc(1, 1, [1, 2, 3, 4, 5; trans_prob_from_prev_state]);
        end
    end
     
end