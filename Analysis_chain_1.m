clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

State_size = 5;
N_chain = 100000;
Time = 100;
pi0 = ones(1,5)/5;
X = chain_1(N_chain, Time, pi0);

%% a) Estimate P(t), and plot the values of its elements over time - for a reasonable range of time.

% Let's estimate transition probability matrix for each value of time.
prob_matrix_estimation = zeros(State_size, State_size, Time - 1);
for time = 1:Time-1
    prob_matrix_estimation(:, :, time) = estimate_transition_matrix(X, time, State_size);
end

% Now, plot how the entries of transition probability matrices change
% over time t.
figure
title('Transition probabilities of chain 1')
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

% Yes, from the plot we can observe that chain 1 is time-homogeneous, as the entries of the transition probability 
% matrix over time t vary only slightly given a large N_chain.

%% c) If it is time-homogeneous:
% c.1) Find its time-homogeneous transition probability and save it as the variable P hat in the 
% file P_hat_chain_i.mat where i is the number of the chain. Choose N_chain and Time in a way to
% have an error (on average) less than 10-3 for each elements of P - explain your approach.


% Since we know that our chain_1 is time-homogenous, thus in order to
% estimate P, we can just set Time=2 with pi0 initalized uniformly, 
% take large N_chain, and then calculate P_hat based on X(1,:) and X(2,:). 
% In order to make sure we have an error (on average) less than 10^-3 for 
% each element of P, we can do this process for a number of times (N times)
% and take P_hat as the average of all estimated P-s.

% Let's take N to be 10.
N=10;

State_size = 5;
% We need very large N_chain to achieve the accuracy we need. Based on
% experiments we saw that 10^7 is more than enough.
N_chain = 10^7;
Time = 2;
pi0 = ones(1,5)/5;

% Let's calculate P_hat N times as explained above.
prob_matrix_estimation = zeros(State_size, State_size, N);
for k = 1:N
    X = chain_1(N_chain, Time, pi0);
    prob_matrix_estimation_cur = estimate_transition_matrix(X, 1, State_size);
    prob_matrix_estimation(:, :, k) = prob_matrix_estimation_cur;
end

accuracy = 10^-3;
size = State_size*State_size;
means = zeros(size,1);
% For making sure that we have an error (on average) less than 10^-3 for
% each element of P, we will use std.
stds = zeros(size,1);

for i = 1:State_size
    for j = 1:State_size
        values = zeros(N, 1);
        for k = 1:N
            values(k) = prob_matrix_estimation(i, j, k);
        end
        means((i-1)*State_size+j) = mean(values);
        stds((i-1)*State_size+j) = std(values);
    end
end

% Let's plot the stds for each entry of P and see that all values are less
% than 10^-3.
figure
title('STD of estimated transtion probabilities for chain 1')
xlabel('Entries of matrix P');
ylabel('Standard Deviation');
hold on;
grid on;
plot(stds,'o');
plot(xlim(),[accuracy,accuracy]);
hold off;

P_hat = transpose(reshape(means,[State_size,State_size]));
% saving P_hat in the file P_hat_chain_1.mat
save('P_hat_chain_1.mat','P_hat');

% loading P_hat from the file P_hat_chain_1.mat
ld = load('P_hat_chain_1.mat');
P_hat = ld.P_hat;

%% c.2) Draw the underlying graph of the chain.

mc = dtmc(P_hat);

figure;
graphplot(mc, 'ColorNodes',true,'ColorEdges',true);

%% d) If it is not time-homogeneous, explain the dynamics of P(t) over time, e.g. how it changes,
% whether it converges, etc.

%% e) When pi0 is an uniform distribution, estimate pi(t) (i.e. the distribution of states at time t,
% i.e. X(t+1,:)), and plot the values of its elements over time.

State_size = 5;
N_chain = 100000;
Time = 100;
pi0 = ones(1,5)/5;
X = chain_1(N_chain, Time, pi0);

% Calculating the distribution of states for each value of time.
pi_t = estimate_distribution(X, Time, State_size);

% Lets draw state distribution over time.
figure
hold on
grid on
title('State distribution of chain 1')
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

% Since our chain is ergodic, i.e. it is irreducible, positive recurrent (follows from having finite
% and irreducible chain) and aperiodic (we have self-loops), hence there
% exists a unique limiting distribution. Thus, we can take any initial
% distribution, say uniform, and find the limiting distribution.

State_size = 5;
N_chain = 10^6;
Time = 100;
pi0 = ones(1,5)/5;
X = chain_1(N_chain, Time, pi0);

% Calculating the distribution of states for each value of time.
pi_t = estimate_distribution(X, Time, State_size);

window_size = 10;
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
title('Max std of state distribution on window of size 10')
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
title('Limiting distribution')
xlabel('State')
ylabel('Probability')   
hold on
grid on
bar(limiting_distr)
hold off

% saving pi_hat in the file pi_hat_chain_1.mat 
save('pi_hat_chain_1.mat','limiting_distr');

% loading pi_hat from the file pi_hat_chain_1.mat
ld = load('pi_hat_chain_1.mat');
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
title('TV for uniform initial dist. for chain 1')
xlabel('Time')
ylabel('Total variation')   
hold on
grid on
set(gca, 'YScale', 'log')
plot(1:Time, total_variation)
hold off

State_size = 5;
N_chain = 10^5;
Time = 100;

pi_ts = zeros(State_size, Time, State_size);
for s = 1:State_size
    pi0 = zeros(1,5);
    pi0(s) = 1;
    X = chain_1(N_chain, Time, pi0);
    pi_ts(:, :, s) = estimate_distribution(X, Time, State_size);
end

total_variations = zeros(Time, State_size);
for s = 1:State_size
    for t = 1:Time
        total_variations(t,s) = sum(abs(pi_ts(:, t, s) - pi_hat)) / 2;
    end
end

figure
title('TV for X(0) = i initial distribution for chain 1')
xlabel('Time')
ylabel('Toral variation')   
hold on
grid on
set(gca, 'YScale', 'log')
colors = ['k','b','r','g','m'];
for s = 1:State_size
    plot(1:Time, total_variations(:,s), 'color', colors(s))
end
legend('X0=1','X0=2','X0=3','X0=4','X0=5');
hold off

% Now, for each of our initial distributions, let's calculate the minimum time needed 
% to have total variation distance less than eps = 0.005.
eps = 0.005;
conv_rates = zeros(State_size,1);
for s = 1:State_size
    for t = 1:Time
        if total_variations(t,s)<eps
            conv_rates(s) = t;
            break;
        end
    end
end

fprintf('X0=%d state has the worst convergence rate\n',find(conv_rates(conv_rates==max(conv_rates))))

T_eps = max(conv_rates);
fprintf('An upper bound for T_eps = %d\n',T_eps);
% An upper bound for T_eps = 62

%% h) If the chain is time-homogeneous, find its stationary distribution using the eigenvalue
% decomposition of your estimation of P. Check whether your findings are consistent with what you 
% found in parts f) and g). Furthermore, plot the stationary distribution with bar-plots.

% Doing eigenvalue decomposition of P_hat. 
% W is a matrix whose columns correspond to the left eigenvectors.
[V, D, W] = eig(P_hat);

% Stationary distribution is the normalized left eigenvector of P_hat
% corresponding to eigenvalue = 1. Hence, we can easily calculate stationary distribution
% from matrix W.
stationary_dstrb = W(:,1) / sum(W(:,1));

figure
title('Stationary distribution')
xlabel('State')
ylabel('Probability')   
hold on
grid on
bar(stationary_dstrb)
hold off

% As we can see, our results are consistent, i.e. our limiting distribution
% found from part f) is the 'same' as the stationary distribution calculated
% from the eigenvalue decomposition of our transition matrix P_hat
fprintf('Euclidean distance between p_hat and estimated stationary distribution = %d\n',...
        sqrt(dot(transpose(pi_hat - stationary_dstrb),(pi_hat - stationary_dstrb))))
% Euclidean distance between p_hat and estimated stationary distribution = 3.942508e-04

% Let's now compute mixing time.
eig_values = diag(D);
alpha_star = max(abs(eig_values(2:State_size)));
spectral_gap = 1 - alpha_star;

mixing_time = - log(2 * eps * min(stationary_dstrb)) / spectral_gap;
% As we can see, our results are consistent with what we found in g)
fprintf('An upper bound for T_eps = %d\n',mixing_time);
% An upper bound for T_eps = 1.045425e+02
