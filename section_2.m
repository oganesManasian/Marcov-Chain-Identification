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

%% a) Write a function as X = MP_chain_i(N_chain, Time, pi_a, x0) which produces N_chain ? N
% realisations of the modified version (in a way to have the limiting distribution pi_a) of chain i with
% length Time ? N and initial state X0 = x0. Note that you are not allowed to directly sample from
% pi_a. The function should be saved as MP_chain_i.m.

% Check applicability of Metropolis-Hasting algorithm ??? Check that matrix
% is symmetric

pi_a = [16,8,4,2,1]/31;
x0 = 1;
Time = 10;
N_chain = 10;
X = MP_chain_1(N_chain, Time, pi_a, x0);

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
Time = 100;
N_chain = 10000;
state_size = 5;

for pi_a_ind = 1:size(pi_a_all, 1)
    display(pi_a_desc(pi_a_ind))
    
    for init_state = 1:5    
        X = MP_chain_1(N_chain, Time, pi_a_all(pi_a_ind, :), init_state);
%         X = MP_chain_1_v2(N_chain, Time, pi_a_all(pi_a_ind, :), init_state);
        
        limiting_distribution = estimate_distribution(X, Time, state_size)';
%         limiting_distribution = estimate_distribution(X, Time, state_size)

        err_mae = mae(pi_a_all(pi_a_ind, :), limiting_distribution);
        fprintf('For initial state %d error is %f\n', init_state, err_mae)
        
    end
end

% Looking at differences between desired distributions and generated ones
% we can see that modeified chain has needed limiting distribution

% TODO Make use of 1F code when implemented for finding limiting distribution

%% b.2) Plot the total variation distance over time, and analyze the effect 
% of the initial state x0 on the convergence rate of the algorithm. Does it
% depend on the desired distribution pi_a? If yes, explain how.

% Compute TV
total_variation = zeros(size(pi_a_all, 1), Time, state_size);

for pi_a_ind = 1:size(pi_a_all, 1)
    for init_state = 1:state_size
        X = MP_chain_1(N_chain, Time, pi_a_all(pi_a_ind, :), init_state);
        
        for time = 1:size(X, 1)
            cur_distribution = estimate_distribution(X, time, state_size)';
            total_variation(pi_a_ind, time, init_state) = ...
                sum(abs(pi_a_all(pi_a_ind, :) - cur_distribution)) / 2; 
        end
    end
end

%% Plot TV
colors = ['k','b','r','g','m'];

figure
for pi_a_ind = 1:size(pi_a_all, 1)
    subplot(3, 1, pi_a_ind)
    title(pi_a_desc(pi_a_ind))
    xlim([0 size(X, 1)])
    xlabel('Time')
    ylabel('Total variation')
    hold on
    grid on
    for init_state = 1:state_size
        plot(1:size(X, 1), total_variation(pi_a_ind, :, init_state), ...
            'DisplayName', sprintf('%i initial state', init_state), ...
            'color', colors(init_state))
    end
end
legend('show');
hold off

% From plot we can see that both desired distribution and initial state
% effect on convergence rate. We can point out two observations: 
% First, the more desired distribution is close to uniform distribution,
% the faster chain converges if we consider average convergence rate among
% all initial states.
% Second, the more state is likely in desired distribution, the faster
% chain converges if we pick this state as initial.

% TODO Make use of 1G code when implemented

%% b.3) For each distribution pi_a, can you estimate (numerically) an upper-bound for Te when
% e = 0.005?
eps = 0.005;

display('Computing mixing time')
for pi_a_ind = 1:size(pi_a_all, 1)
    display(pi_a_desc(pi_a_ind))
    
    for time = 1:size(total_variation, 2)
        if max(total_variation(pi_a_ind, time, :)) < eps
            display(sprintf('Mixing time is %d', time))
            break
        end
        
        if time == size(total_variation, 2)
            display('Mixing time was not found. Please increase considering time')
            fprintf('TV on max time is %d\n', max(total_variation(pi_a_ind, time, :)))
        end
    end
end

%% c) For each of the three choices of pi_a in part b), which chain has a better convergence rate? Can 
% you explain your observations intuitively?



