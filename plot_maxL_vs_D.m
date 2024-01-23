%% plot MaxL leakage vs D in different p

clc
clear all 
close all

% assign value
N = 3;
K = 2;
r_diff = [0.1 0.3 0.6];
r_same = [0.333 0.333 0.333];

% calculate D*
D_star = 0;
for i = 1:K
    D_star = D_star + (1/N)^(i-1);
end

D = 1.00:0.02:D_star;

% create initial matirx 
MaxL_leakage_1 = zeros(size(D));
MaxL_leakage_2 = zeros(size(D));
MaxL_leakage_3 = zeros(size(D));
MaxL_leakage_4 = zeros(size(D));
MaxL_leakage_5 = zeros(size(D));

% %% homo
% for d = 1:length(D)
%     % line 1 : weight same optimum without #
%     [MaxL_leakage_1(d), p__rsame_opt] = auto_compute_maxL_adjustp(N,K,r_same,D(d),0,0);
% 
%     % line 2 : weight same optimum with # and symmetric allocation
%     [MaxL_leakage_2(d), p__rsame_opt] = auto_compute_maxL_adjustp(N,K,r_same,D(d),5,1);
% 
%     % line 3 : weight same optimum with # and allocation we found
%     [MaxL_leakage_3(d), p__rsame_opt] = auto_compute_maxL_adjustp(N,K,r_same,D(d),4,1);
% end
% 
% %plot
% figure(1)
% plot(MaxL_leakage_1,D)
% hold on
% plot(MaxL_leakage_2,D,'-o')
% hold on
% plot(MaxL_leakage_3,D)
% % hold on
% % plot(MaxL_leakage_4,D)
% % hold on
% % plot(MaxL_leakage_5,D)
% xlabel('$\rho$','Interpreter','latex')
% ylabel('D')
% legend({"opt w/o p# ", "symmetric allocation w/ p#", "paper allocation w/ p#"})
% title("N = 3, K = 2 Homogeneous case ")

%% hetero
for d = 1:length(D)
    % line 1 : weight same optimum without #
    [MaxL_leakage_1(d), p__rsame_opt] = auto_compute_maxL_adjustp(N,K,r_diff,D(d),0,0);
    
    % line 2 : weight same optimum with # and symmetric allocation
    [MaxL_leakage_2(d), p__rsame_opt] = auto_compute_maxL_adjustp(N,K,r_diff,D(d),5,1);

    % line 3 : weight same optimum with # and allocation we found
    [MaxL_leakage_3(d), p__rsame_opt] = auto_compute_maxL_adjustp(N,K,r_diff,D(d),4,1);
end

%plot
figure(1)
plot(MaxL_leakage_1,D)
hold on
plot(MaxL_leakage_2,D,'-o')
hold on
plot(MaxL_leakage_3,D)
% hold on
% plot(MaxL_leakage_4,D)
% hold on
% plot(MaxL_leakage_5,D)
xlabel('$\rho$','Interpreter','latex')
ylabel('D')
legend({"opt w/o p# ", "symmetric allocation w/ p#", "paper allocation w/ p#"})
title("N = 3, K = 2 Heterogeneous case ")