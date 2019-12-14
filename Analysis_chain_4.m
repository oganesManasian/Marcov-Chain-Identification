clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

State_size = 5;
N_chain = 100000;
Time = 200;
pi0 = ones(1,5)/5;
X = chain_4(N_chain, Time, pi0);
%% a) Estimate P(t), and plot the values of its elements over time - for a reasonable range of time.

% Let's estimate transition probability matrix for each value of time.
prob_matrix_estimation = zeros(State_size, State_size, Time - 1);
for time = 1:Time-1
    prob_matrix_estimation(:, :, time) = estimate_transition_matrix(X, time, State_size);
end

% Now, plot how the entries of transition probability matrices change
% over time t.
figure
title('Transition probabilities of chain 4')
xlabel('Time')
ylabel('Probability')   
hold on
grid on

for i = 1:State_size
    for j = 1:State_size
        plot(1:Time-1,  squeeze(prob_matrix_estimation(i, j, :)))
    end
end

hold off
%% b) Is this chain time-homogeneous?

% As the transition matrix changes in time, the chain is not
% time-homogenious.

%% d) If it is not time-homogeneous, explain the dynamics of P(t) over time, e.g. how it changes,
% whether it converges, etc.

% From the plot we see that each element of the transition matrix converges
% to some value. Let's plot the total variation from the final matrix in time.

limiting_transition_matrix = mean(prob_matrix_estimation(:, :, Time-5:Time-1), 3);

total_variation = zeros(Time - 1);
for time = 1:Time-1
    total_variation(time) = sum(abs(limiting_transition_matrix - ... 
                                    prob_matrix_estimation(:, :, time)), 'ALL') / 2;
end

figure
title('TV of transition matrix and limiting matrix chain 4')
xlabel('Time')
ylabel('Probability')   
hold on
grid on

plot(1:Time-1,  total_variation);

hold off

mc = dtmc(limiting_transition_matrix);
figure;
graphplot(mc, 'ColorNodes',true,'ColorEdges',true);

% Indeed, we see that the total distance is exponentially decreasing in time!
%% e) When pi0 is an uniform distribution, estimate pi(t) (i.e. the distribution of states at time t,
% i.e. X(t+1,:)), and plot the values of its elements over time.

% Calculating the distribution of states for each value of time.
pi_t = estimate_distribution(X, Time, State_size);

% Lets draw state distribution over time.
figure
hold on
grid on
title('State distribution of chain 4')
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

% Even though the chain is not time-homogenious, we know that its
% transition matrix converges. So, to say if there is a limiting distribution 
% or not we can study the limiting matrix as after some time 
% the chain  with any needed precision is time-homogeniuos. 
% Seeing the limiting matrix we can notice that it is ergodic, hence there
% exists a unique limiting distribution. Thus, we can take any initial
% distribution, say uniform, and find the limiting distribution.

N_chain = 10^6;
Time = 200;
X = chain_4(N_chain, Time, pi0);

% Calculating the distribution of states for each value of time.
pi_t = estimate_distribution(X, Time, State_size);

window_size = 20;
threshold = 10^-3;

errors = zeros(Time - window_size, 1);
state_errors = zeros(State_size, Time - window_size);

for t = 1:Time - window_size
    for state = 1:State_size
        state_errors(state, t) = std(pi_t(state, t:t + window_size));
    end
    % taking the maximum of the variation over our states
    errors(t) = max(state_errors(:, t));
end

figure
title('Max std of state distribution on 10-size window for chain 4')
xlabel('Time')
ylabel('STD')   
hold on
grid on
set(gca, 'YScale', 'log')
plot(1+window_size:Time, errors)
hold off

limiting_distr = zeros(State_size, 1);
limiting_distr_slice = zeros(State_size, window_size);
for t = 1:Time - window_size
    if errors(t) < threshold
        limiting_distr_slice = pi_t(:, t:t + window_size);
        limiting_distr = mean(limiting_distr_slice, 2); 
    end
end

figure
title('Limiting distribution of chain 4')
xlabel('State')
ylabel('Probability')   
hold on
grid on
bar(limiting_distr)
hold off

% saving pi_hat in the file pi_hat_chain_4.mat 
save('pi_hat_chain_4.mat','limiting_distr');

% loading pi_hat from the file pi_hat_chain_4.mat
ld = load('pi_hat_chain_4.mat');
pi_hat = ld.limiting_distr;
%% g) For the case where the limiting distribution pi exists, plot the total-variation distance over
% time. Plot the total-variation distance over time also for the cases where X(0) = i, 
% for each i in {1, 2, 3, 4, 5}. Which initial state has the worst convergence rate? Can you estimate 
% (numerically) an upperbound for Te when e = 0.005?

% Let's calculate total-variation distance for each value of time. 
total_variation = zeros(Time, 1);
for t = 1:Time
    total_variation(t) = sum(abs(pi_t(:, t) - pi_hat)) / 2;
end

figure
title('TV for uniform initial dist. for chain 4')
xlabel('Time')
ylabel('Total variation')   
hold on
grid on
set(gca, 'YScale', 'log')
plot(1:Time, total_variation)
hold off

% And for different initial distributions

State_size = 5;
N_chain = 10^5;
Time = 200;

pi_ts = zeros(State_size, Time, State_size);
for s = 1:State_size
    pi0 = zeros(1,5);
    pi0(s) = 1;
    X_tmp = chain_4(N_chain, Time, pi0);
    pi_ts(:, :, s) = estimate_distribution(X_tmp, Time, State_size);
end

total_variation_for_states = zeros(Time, State_size);
for s = 1:State_size
    for t = 1:Time
        total_variation_for_states(t,s) = sum(abs(pi_ts(:, t, s) - pi_hat)) / 2;
    end
end

figure
title('TV for X(0) = i initial distribution for chain 4')
xlabel('Time')
ylabel('Total variation')   
hold on
grid on
set(gca, 'YScale', 'log')
colors = ['k','b','r','g','m'];
for s = 1:State_size
    plot(1:Time, total_variation_for_states(:,s), 'color', colors(s))
end
legend('X0=1','X0=2','X0=3','X0=4','X0=5');
hold off

% Now, for each of our initial distributions, let's calculate the minimum time needed 
% to have total variation distance less than eps = 0.005.
eps = 0.005;
conv_rates = zeros(State_size,1);
for s = 1:State_size
    for t = 1:Time
        if total_variation_for_states(t,s)<eps
            conv_rates(s) = t;
            break;
        end
    end
end

fprintf('X0=%d state has the worst convergence rate\n',find(conv_rates(conv_rates==max(conv_rates))))
% X0=1 state has the worst convergence rate
T_eps = max(conv_rates);
fprintf('An upper bound for T_eps = %d\n',T_eps);
% An upper bound for T_eps = 105