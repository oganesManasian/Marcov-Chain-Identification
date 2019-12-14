clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

% Assume that chain i is time-homogeneous. The aim of this part is to use the Metropolis-Hasting
% algorithm to change this chain in a way to have the arbitrary distribution pi_a as its limiting distribution.
% Given Xt as the state of the chain at time t, the next state Xt+1 can be sampled as [Xt; Xt+1]
% = chain i(1, 1, Xt). With the Metropolis-Hasting algorithm, you can find the probability at+1
% with which you should accept the new sample Xt+1 (otherwise stay at the same state as Xt) to
% finally have a limiting distribution pi_a.
% For each time-homogeneous Markov Chain among the 4 chains (i.e. chain 1.p, chain 2.p, chain 3.p,
% and chain 4.p)

% Chains 1, 2 are time-homogeneous, let's analyze them
% Let's analyse chain 2

%% a) Write a function as X = MP_chain_i(N_chain, Time, pi_a, x0) which produces N_chain ? N
% realisations of the modified version (in a way to have the limiting distribution pi_a) of chain i with
% length Time ? N and initial state X0 = x0. Note that you are not allowed to directly sample from
% pi_a. The function should be saved as MP_chain_i.m.

%%  b) For each of your chains, and for each of the cases
% 1. pi_a = [16,8,4,2,1]/31;
% 2. pi_a = [1,1,4,1,1]/8;
% 3. pi_a = [4,2,1,2,4]/13;

pi_a_desc = {'Desired distribution: [16,8,4,2,1]/31', ...
             'Desired distribution: [1,1,4,1,1]/8', ...
             'Desired distribution: [4,2,1,2,4]/13'};
pi_a_all = [[16,8,4,2,1]/31; [1,1,4,1,1]/8; [4,2,1,2,4]/13];

% do the following steps:

%% b.1) Evaluate your code, i.e. show that your modified chain has a limiting distribution equal to
% pi_a for all choices of the initial state x0.

Time = 150;
N_chain = 10^5;
state_size = 5;

% Compute TV
total_variation = zeros(size(pi_a_all, 1), Time, state_size);

for pi_a_ind = 1:size(pi_a_all, 1)
    disp(pi_a_desc(pi_a_ind))
    pi_a = pi_a_all(pi_a_ind, :);
    
    for init_state = 1:state_size
        X = MP_chain_2(N_chain, Time, pi_a, init_state);
        
        estimated_distrubution = estimate_distribution(X, Time, state_size);
        total_variation(pi_a_ind, :, init_state) = ...
                sum(abs(repmat(pi_a', 1, Time) - estimated_distrubution)) / 2; 
        
        error = total_variation(pi_a_ind, Time, init_state);
        fprintf('For initial state %d error is %f\n', init_state, error)
    end
end

% OUTPUT:
% 'Desired distribution: [16,8,4,2,1]/31'
% For initial state 1 error is 0.002423
% For initial state 2 error is 0.002483
% For initial state 3 error is 0.001677
% For initial state 4 error is 0.001931
% For initial state 5 error is 0.002079
%     'Desired distribution: [1,1,4,1,1]/8'
% For initial state 1 error is 0.002040
% For initial state 2 error is 0.003260
% For initial state 3 error is 0.005190
% For initial state 4 error is 0.002390
% For initial state 5 error is 0.001660
%     'Desired distribution: [4,2,1,2,4]/13'
% For initial state 1 error is 0.000952
% For initial state 2 error is 0.002972
% For initial state 3 error is 0.002562
% For initial state 4 error is 0.002162
% For initial state 5 error is 0.003775

% Value of total variation at time 150 are very small for each
% initial state and thus pi_a are indeed limiting distributions

%% b.2) Plot the total variation distance over time, and analyze the effect 
% of the initial state x0 on the convergence rate of the algorithm. Does it
% depend on the desired distribution pi_a? If yes, explain how.

%% Plot TV
colors = ['k','b','r','g','m'];

figure
for pi_a_ind = 1:size(pi_a_all, 1)
    subplot(3, 1, pi_a_ind) 
    title(pi_a_desc(pi_a_ind))
    xlim([1 min(60, size(X, 1))]) % After Time=60 curve is very close to 0
    xlabel('Time')
    ylabel('Total variation')
    set(gca, 'YScale', 'log')
    hold on
    grid on
    for init_state = 1:state_size
        plot(1:size(X, 1), total_variation(pi_a_ind, :, init_state), ...
            'displayName', sprintf('%i initial state', init_state), ...
            'color', colors(init_state))
    end
end
legend('show');
hold off

% From the plot we can see that both desired distribution and initial state
% affect the convergence rate. There are two main observations: 
% First, since base chain has non-zero transition probabilities only for
% neighboring states, new chain converges faster when in desired
% distribution most probable states are neighboring. As an instances we see
% how fast the chain converges in option 1 of desired distribution where most
% probable states are 1 and 2. Conversely, mixing time of chain with option
% 3 of desired distribution, where most probable states are 1 and 5,
% is the highest.
% Second, the higher probability of the state in desired distribution, the faster
% chain converges if we pick this state as initial.

%% Extra plot. Total variation when uniform distribution is initial
Time = 150;
N_chain = 10^5;
state_size = 5;

x0 = randsrc(1, N_chain, [1, 2, 3, 4, 5; 0.2, 0.2, 0.2, 0.2, 0.2]);

total_variation_uniform = zeros(size(pi_a_all, 1), Time);

for pi_a_ind = 1:size(pi_a_all, 1)
    disp(pi_a_desc(pi_a_ind))
    pi_a = pi_a_all(pi_a_ind, :);
    
    X = MP_chain_2(N_chain, Time, pi_a, x0);
        
    estimated_distrubution = estimate_distribution(X, Time, state_size);
    total_variation_uniform(pi_a_ind, :) = ...
        sum(abs(repmat(pi_a', 1, Time) - estimated_distrubution)) / 2; 
end

%% Extra plot. Total variation when dominant distribution is initial
dominant_value = 0.9;
recessive_value = (1 - dominant_value) / (state_size - 1);
total_variation_dominant = zeros(size(pi_a_all, 1), Time, state_size);

for pi_a_ind = 1:size(pi_a_all, 1)
    disp(pi_a_desc(pi_a_ind))
    pi_a = pi_a_all(pi_a_ind, :);
    
    for dominant_state = 1:state_size
        distr_matr = [1, 2, 3, 4, 5; ...
            recessive_value, recessive_value, recessive_value, recessive_value, recessive_value];
        distr_matr(2, dominant_state) = dominant_value;
        x0 = randsrc(1, N_chain, distr_matr);
    
        X = MP_chain_2(N_chain, Time, pi_a, x0);
        
        estimated_distrubution = estimate_distribution(X, Time, state_size);
        total_variation_dominant(pi_a_ind, :, dominant_state) = ...
                sum(abs(repmat(pi_a', 1, Time) - estimated_distrubution)) / 2; 
    end
end

%% Extra plot 1.
colors = ['k','b','r','g','m'];

figure
title('Extra plot. Total Variation for different initial distributions')
for pi_a_ind = 1:size(pi_a_all, 1)
    subplot(3, 1, pi_a_ind) 
    title(pi_a_desc(pi_a_ind))
    xlim([1 min(60, size(X, 1))]) % After Time=60 curve is very close to 0
    xlabel('Time')
    ylabel('Total variation')
%     set(gca, 'YScale', 'log')
    hold on
    grid on

    plot(1:size(X, 1), total_variation_uniform(pi_a_ind, :), ...
        'displayName', 'Uniform', 'color', 'k', 'lineStyle', '^')
    for dominant_state = 1:state_size
        plot(1:size(X, 1), total_variation_dominant(pi_a_ind, :, dominant_state), ...
            'displayName', sprintf('%i as dominant state', dominant_state), ...
            'color', colors(dominant_state))
    end
end
legend('show');
hold off

%% Extra plot 2.
colors = ['k','b','r','g','m'];

figure
title('Extra plot. Total Variation when uniform is initial distribution')
xlim([1 min(60, size(X, 1))]) % After Time=60 curve is very close to 0
xlabel('Time')
ylabel('Total variation')
% set(gca, 'YScale', 'log')
hold on
grid on
    
for pi_a_ind = 1:size(pi_a_all, 1)
    plot(1:size(X, 1), total_variation_uniform(pi_a_ind, :), ...
        'displayName', pi_a_desc(pi_a_ind), ...
        'color', colors(pi_a_ind), 'lineStyle', '^')
end
legend('show');
hold off

%% b.3) For each distribution pi_a, can you estimate (numerically) an upper-bound for Te when
% e = 0.005?
eps = 0.005;

fprintf('Computing mixing time\n')
for pi_a_ind = 1:size(pi_a_all, 1)
    disp(pi_a_desc(pi_a_ind))
    
    for time = 1:size(total_variation, 2)
        if max(total_variation(pi_a_ind, time, :)) < eps
            fprintf(sprintf('Mixing time is %d\n', time))
            break
        end
        
        if time == size(total_variation, 2)
            fprintf('Mixing time was not found. Please increase considering time')
            fprintf('TV on max time is %d\n', max(total_variation(pi_a_ind, time, :)))
        end
    end
end

% 'Desired distribution: [16,8,4,2,1]/31' Mixing time is 31    
% 'Desired distribution: [1,1,4,1,1]/8'   Mixing time is 27    
% 'Desired distribution: [4,2,1,2,4]/13'  Mixing time is 93

%% c) For each of the three choices of pi_a in part b), which chain has a 
% better convergence rate? Can you explain your observations intuitively?



