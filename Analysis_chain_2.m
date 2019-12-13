clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

State_size = 5;
N_chain = 100000;
Time = 100;
pi0 = ones(1,5)/5;
X = chain_2(N_chain, Time, pi0);

%% a) Estimate P(t), and plot the values of its elements over time - for a reasonable range of time.

% Let's estimate transition probability matrix for each value of time.
prob_matrix_estimation = zeros(State_size, State_size, Time - 1);
for time = 1:Time-1
    prob_matrix_estimation(:, :, time) = estimate_transition_matrix(X, time, State_size);
end

% Now, plot how the entries of transition probability matrices change
% over time t.
figure
title('Transition probabilities of chain 2')
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

% Yes, from the plot we can observe that chain 2 is time-homogeneous, as the entries of the transition probability 
% matrix over time t vary only slightly given a large N_chain.

%% c) If it is time-homogeneous:
% c.1) Find its time-homogeneous transition probability and save it as the variable P hat in the 
% file P_hat_chain_i.mat where i is the number of the chain. Choose N_chain and Time in a way to
% have an error (on average) less than 10-3 for each elements of P - explain your approach.


% Since we know that our chain_2 is time-homogenous, thus in order to
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
    X = chain_2(N_chain, Time, pi0);
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
xlabel('Entries of matrix P for chain 2');
ylabel('Standard Deviation');
hold on;
grid on;
plot(stds,'o');
plot(xlim(),[accuracy,accuracy]);
hold off;

P_hat = transpose(reshape(means,[State_size,State_size]));
% saving P_hat in the file P_hat_chain_2.mat
save('P_hat_chain_2.mat','P_hat');

% loading P_hat from the file P_hat_chain_2.mat
ld = load('P_hat_chain_2.mat');
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
X = chain_2(N_chain, Time, pi0);

% Calculating the distribution of states for each value of time.
pi_t = estimate_distribution(X, Time, State_size);

% Lets draw state distribution over time.
figure
hold on
grid on
title('State distribution of chain 2')
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

% As we can observe from the plot done in e), there is no limiting
% distribution. After some time t0, if look at the sequence of the state distributions t>t0, 
% we can see that for even values of t and for odd values of t, the
% limit distribution is different, meaning that our sequence has no limit,
% i.e. no limiting distribution exists.
%% g) For the case where the limiting distribution pi exists, plot the total-variation distance over
% time. Plot the total-variation distance over time also for the cases where X(0) = i, 
% for each i in {1, 2, 3, 4, 5}. Which initial state has the worst convergence rate? Can you estimate 
% (numerically) an upperbound for Te when e = 0.005?

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
title('Stationary distribution of chain 2')
xlabel('State')
ylabel('Probability')   
hold on
grid on
bar(stationary_dstrb)
hold off
