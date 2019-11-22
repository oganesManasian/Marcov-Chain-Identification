clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

State_size = 5;
N_chain = 10000;
Time = 100;
pi0 = ones(1,5)/5;
X = chain_1(N_chain, Time, pi0);

%% a) Estimate P(t), and plot the values of its elements over time - for a reasonable range of time.

% Lets try to estimate probability matrix for each value of time
prob_matrix_estimation = zeros(State_size, State_size, Time - 1);
for time = 1:Time-1
%     fprintf('Time %i', time)
    prob_matrix_estimation_cur = estimate_transition_matrix(X, time, State_size);
%     display(prob_matrix_estimation_cur)
    prob_matrix_estimation(:, :, time) = prob_matrix_estimation_cur;
end

% Now we have to plot changes of probability of transition for each path
% over time

figure
title('Transition probabilities')
xlabel('Time')
ylabel('Probability')   
hold on
grid on

for i = 1:State_size
    for j = 1:State_size
        values = zeros(Time - 1, 1);
        for k = 1:Time - 1
            values(k) = prob_matrix_estimation(i, j, k);
        end
%         plot(1:Time-1, values, 'DisplayName', sprintf('%i->%i', i, j))
        plot(1:Time-1, values)
    end
end

% legend('show')
hold off

%% b) Is this chain time-homogeneous?

% No, plot reveals that values of transition probabilities matrix changes
% over time

%% c) If it is time-homogeneous:
% c.1) Find its time-homogeneous transition probability and save it as the variable P hat in the 
% file P_hat_chain_i.mat where i is the number of the chain. Choose N_chain and Time in a way to
% have an error (on average) less than 10?3 for each elements of P - explain your approach.
% c.2) Draw the underlying graph of the chain.

% It's not time-homogeneous

%% d) If it is not time-homogeneous, explain the dynamics of P(t) over time, e.g. how it changes,
% whether it converges, etc.

% Yes, they converge (some of them after time = 40, some of them much
% earlier). TODO: rewrite

%% e) When pi0 is an uniform distribution, estimate pi(t) (i.e. the distribution of states at time t,
% i.e. X(t+1,:)), and plot the values of its elements over time.

pi_t = zeros(State_size, Time);
for time = 1:Time
%     fprintf('Time %i', time)
    pi_t(:, time) = estimate_distribution(X, time, State_size);
%     display(pi_t)
end

% Lets draw state distribution over time
figure
hold on
grid on
title('State distribution')
xlabel('Time')
ylabel('Probability')
colors = ['k','b','r','g','m'];
for state = 1:State_size
    plot(1:Time, pi_t(state, :), 'DisplayName', sprintf('%i state', state), 'color', colors(state))
end

legend('show');
hold off

%% f) Does there exist any limiting distribution for this chain? If yes, save it as the variable 
% pi_hat in the file pi_hat_chain_i.mat where i is the number of the chain. Choose N_chain and Time 
% in a way to have an error (on average) less than 10-3 for each element of pi - explain your 
% approach. Furthermore, plot the values of the limiting distribution as a bar-plot.

%% g) For the case where the limiting distribution pi exists, plot the total-variation distance over
% time. Plot the total-variation distance over time also for the cases where X(0) = i, 
% for each i in {1, 2, 3, 4, 5}. Which initial state has the worst convergence rate? Can you estimate 
% (numercically) an upperbound for Te when e = 0.005?

%% h) If the chain is time-homogeneous, find its stationary distribution using the eigenvalue
% decomposition of your estimation of P. Check whether your findings are consistent with what you 
% found in parts f) and g). Furthermore, plot the stationary distribution with bar-plots.
