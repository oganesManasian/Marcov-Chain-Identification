clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
clear;  % Erase all existing variables.
workspace;  % Make sure the workspace panel is showing.

State_size = 5;
N_chain = 100000;
Time = 100;
pi0 = ones(1,5)/5;
X = chain_3(N_chain, Time, pi0);
%% a) Estimate P(t), and plot the values of its elements over time - for a reasonable range of time.

% Let's estimate transition probability matrix for each value of time.
prob_matrix_estimation = zeros(State_size, State_size, Time - 1);
for time = 1:Time-1
    prob_matrix_estimation(:, :, time) = estimate_transition_matrix(X, time, State_size);
end

% Now, plot how the entries of transition probability matrices change
% over time t.
figure
title('Transition probabilities of chain 3')
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

% No, we can clearly see that the chain is not time-homogeneous.
% Moreover its transitions matrix seems to be periodic
% and we can sugest that the period is 2, but we will check it later

%% d) If it is not time-homogeneous, explain the dynamics of P(t) over time, e.g. how it changes,
% whether it converges, etc.

% As we suspect from the plot, transition matrix is 2-periodic. (by periodicyty per)
% To check this we can plot the transition matrix for odd and even steps.


figure
title('Transition probabilities of chain 3 for odd times')
xlabel('Time')
ylabel('Probability')   
hold on
grid on

for i = 1:State_size
    for j = 1:State_size
        plot(1:2:Time-1, squeeze(prob_matrix_estimation(i, j, 1:2:Time-1)))
    end
end

hold off


figure
title('Transition probabilities of chain 3 for even times')
xlabel('Time')
ylabel('Probability')   
hold on
grid on

for i = 1:State_size
    for j = 1:State_size
        plot(2:2:Time-1, squeeze(prob_matrix_estimation(i, j, 2:2:Time-1)))
    end
end

hold off

% Indeed, our intuition was right – the transition matrix 
% is alternation of twomatrices P1, P2, P1, P2 ...
%% e) When pi0 is an uniform distribution, estimate pi(t) (i.e. the distribution of states at time t,
% i.e. X(t+1,:)), and plot the values of its elements over time.

% Calculating the distribution of states for each value of time.
pi_t = estimate_distribution(X, Time, State_size);

% Lets draw state distribution over time.
figure
hold on
grid on
title('State distribution in time for chain 3')
xlabel('Time')
ylabel('Probability')
colors = ['k','b','r','g','m'];
for state = 1:State_size
    plot(1:Time, pi_t(state, :), 'DisplayName', sprintf('%i state', state), 'color', colors(state))
end

legend('show');
hold off
%% f) Does there exist any limiting distribution for this chain? If yes, save it as the variable 
% We can clearly see from the plot – there is no limiting destribution
