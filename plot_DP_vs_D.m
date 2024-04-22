%% plot Download cost vs DP constraint in different epsilons

clc
clear all 
close all

% assign value
N = 20; % number of servers
K = 5; % number of messages
r_diff = [0.1 0.3 0.6]; % homo weights
r_same = [0.333 0.333 0.333]; % hetero weights

% Choice of different epsilons
num_epoch = 200; % number of epochs
eps = linspace(0, 20, num_epoch);

% calculate lower/upper bound of Download cost
D_low = zeros(size(eps)); % lower bound of D
D_up = zeros(size(eps)); % upper bound of D
for d = 1:length(eps)
    for i = 1:K
        D_low(d) = D_low(d) + (N*exp(eps(d)))^(1-i);
    end
    D_up(d) = 1 + (N^(K-1)-1) / ((N-1)*(exp(eps(d))+N^(K-1)-1));
end

% Create matrix to store D
DP_cost_homo = zeros(size(eps)); % Download cost D

%% Numerically solve optimal D vs DP leakage
for d = 1:length(eps)
    % [DP_cost_homo(d)] = auto_compute_minD_DP(N, K,eps(d));
    DP_cost_homo(d) = N/(N-1) - exp(eps(d)*(K-1)) / (N-1) / (exp(eps(d))+N-1)^(K-1);
    % [DP_cost_homo(d)] = auto_compute_minD_DP_LPreduced(N, K, eps(d));
end

%% Calculate the gap
Gap_TSC_low = DP_cost_homo ./ D_low;
Gap_Ravi_low = D_up ./ D_low;
Gap_Ravi_TSC = D_up ./ DP_cost_homo;

%% Plot the result
figure
set(gcf,'position',[50, 100, 550, 400])
set(gcf,'Color',[0.9,0.9,0.9])
plot(eps,D_up,'LineWidth',2)
hold on
plot(eps,D_low,'LineWidth',2)
hold on
plot(eps,DP_cost_homo,'LineWidth',2)
hold on
plot(eps,Gap_TSC_low,'--','LineWidth',1.5)
hold on
% plot(eps,Gap_Ravi_low,'--','LineWidth',1.5)
% hold on
plot(eps,Gap_Ravi_TSC,'--','LineWidth',1.5)

grid on
xlabel('User Privacy Leakage ($\epsilon$)','Interpreter','latex')
% ylim([1 1.021])
ylabel('Download Cost $D$','Interpreter','latex')
legend("$D^{UB}(\epsilon)$", "$D^{LB}(\epsilon)$", "$D(\epsilon)$", "$D(\epsilon) / D^{LB}(\epsilon)$", "$D^{UB}(\epsilon) / D(\epsilon)$", 'Interpreter','latex')
title("N = 50, K = 50")