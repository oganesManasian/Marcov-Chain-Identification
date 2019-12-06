% Assume that chain i is time-homogeneous. The aim of this part is to use the Metropolis-Hasting
% algorithm to change this chain in a way to have the arbitrary distribution pi_a as its limiting distribution.
% Given Xt as the state of the chain at time t, the next state Xt+1 can be sampled as [Xt; Xt+1]
% = chain i(1, 1, Xt). With the Metropolis-Hasting algorithm, you can find the probability at+1
% with which you should accept the new sample Xt+1 (otherwise stay at the same state as Xt) to
% finally have a limiting distribution pi_a.
% For each time-homogeneous Markov Chain among the 4 chains (i.e. chain 1.p, chain 2.p, chain 3.p,
% and chain 4.p)

% Chains 1,... are time-homogeneous
Xprev = 1;
pair = chain_1(1, 2, Xprev);
Xnew = pair(2);
% Finding acceptance rate




%% a) Write a function as X = MP_chain_i(N_chain, Time, pi_a, x0) which produces N_chain ? N
% realisations of the modified version (in a way to have the limiting distribution pi_a) of chain i with
% length Time ? N and initial state X0 = x0. Note that you are not allowed to directly sample from
% pi_a. The function should be saved as MP_chain_i.m.

% Check applicability of Metropolis-Hasting algorithm ??? Check that matrix
% is symmetric

pi_a = [16,8,4,2,1]/31;
x0 = 1;
X = MP_chain_1(10, 10, pi_a, x0);

%%  b) For each of your chains, and for each of the cases
% 1. pi_a = [16,8,4,2,1]/31;
% 2. pi_a = [1,1,4,1,1]/8;
% 3. pi_a = [4,2,1,2,4]/13;
% do the following steps:



%% b.1) Evaluate your code, i.e. show that your modified chain has a limiting distribution equal to
% pi_a for all choices of the initial state x0.


% Make use of 1G code when implemented

%% b.2) Plot the total variation distance over time, and analyze the effect of the initial state x0 on the
% convergence rate of the algorithm. Does it depend on the desired distribution pi_a? If yes, explain
% how.



%% b.3) For each distribution pi_a, can you estimate (numerically) an upper-bound for Te when
% e = 0.005?


%% c) For each of the three choices of pi_a in p


