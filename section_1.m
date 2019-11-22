N_chain = 5;
Time = 2;
pi0 = ones(1,5)/5;
X = chain_1(N_chain, Time, pi0);


% a) Estimate P(t), and plot the values of its elements over time - for a reasonable range of time.

% b) Is this chain time-homogeneous?

% c) If it is time-homogeneous:
% c.1) Find its time-homogeneous transition probability and save it as the variable P hat in the 
% file P_hat_chain_i.mat where i is the number of the chain. Choose N_chain and Time in a way to
% have an error (on average) less than 10?3 for each elements of P - explain your approach.
% c.2) Draw the underlying graph of the chain.

% d) If it is not time-homogeneous, explain the dynamics of P(t) over time, e.g. how it changes,
% whether it converges, etc.

% e) When pi0 is an uniform distribution, estimate ?(t) (i.e. the distribution of states at time t,
% i.e. X(t+1,:)), and plot the values of its elements over time.

% f) Does there exist any limiting distribution for this chain? If yes, save it as the variable 
% pi_hat in the file pi_hat_chain_i.mat where i is the number of the chain. Choose N_chain and Time 
% in a way to have an error (on average) less than 10?3 for each element of ? - explain your 
% approach. Furthermore, plot the values of the limiting distribution as a bar-plot.

% g) For the case where the limiting distribution ? exists, plot the total-variation distance over
% time. Plot the total-variation distance over time also for the cases where X(0) = i, 
% ?i ? {1, 2, 3, 4, 5}. Which initial state has the worst convergence rate? Can you estimate 
% (numercically) an upperbound for Te when e = 0.005?

% h) If the chain is time-homogeneous, find its stationary distribution using the eigenvalue
% decomposition of your estimation of P. Check whether your findings are consistent with what you 
% found in parts f) and g). Furthermore, plot the stationary distribution with bar-plots.
