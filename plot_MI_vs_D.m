%% plot MI leakage vs D in different p

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
MI_leakage_homo_1 = zeros(size(D));
MI_leakage_homo_2 = zeros(size(D));
MI_leakage_homo_3 = zeros(size(D));
MI_leakage_homo_4 = zeros(size(D));
MI_leakage_homo_5 = zeros(size(D));

MI_leakage_hete_1 = zeros(size(D));
MI_leakage_hete_2 = zeros(size(D));
MI_leakage_hete_3 = zeros(size(D));
MI_leakage_hete_4 = zeros(size(D));
MI_leakage_hete_5 = zeros(size(D));

%% Homo
for d = 1:length(D)
    % line 1 : optimal
    [MI_leakage_homo_1(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_same,D(d),0,1);

    % line 2 : no p#
    [MI_leakage_homo_2(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_same,D(d),0,0);

end

%% Hetero
for d = 1:length(D)
    % line 1 : optimal
    [MI_leakage_hete_1(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),0,1);

    % line 2 : no p#
    [MI_leakage_hete_2(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),0,0);
    
    % line 3 :pw uniform
    [MI_leakage_hete_3(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),3,1);

    % % line 4 :pw uniform, p# symmetric
    % [MI_leakage_4(d), p__rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),8,1);
end

%% plot homo
figure
set(gcf,'position',[50, 100, 550, 400])
set(gcf,'Color',[0.9,0.9,0.9])
plot(MI_leakage_homo_2,D,'-square','LineWidth',1.2)
hold on
plot(MI_leakage_homo_1,D,'-o','LineWidth',1.2)
grid on
xlabel('MI constraint $\rho$','Interpreter','latex')
ylabel('Download cost $D$','Interpreter','latex')
legend("Optimal allocation w/o p# ", "Paper allocation w/ p#")
title("Homogeneous case: N = 3, K = 2")

%% plot hetero
% figure
% plot(MI_leakage_1,D,'-x')
% hold on
% plot(MI_leakage_2,D,'-o')
% hold on
% plot(MI_leakage_3,D)
% % hold on
% % plot(MaxL_leakage_4,D)
% xlabel('$\rho$','Interpreter','latex')
% ylabel('D')
% legend({"opt w/ p# ", "opt w/o p#","pw uniform"})
% title("N = 3, K = 2 Heterogeneous case ")


figure
set(gcf,'position',[50, 100, 550, 400])
set(gcf,'Color',[0.9,0.9,0.9])
plot(MI_leakage_hete_2,D,'-square','LineWidth',1.2)
hold on
plot(MI_leakage_hete_3,D,'-o','LineWidth',1.2)
hold on
plot(MI_leakage_hete_1,D,'-*','LineWidth',1.2)
grid on
xlabel('MI constraint $\rho$','Interpreter','latex')
ylabel('Download cost $D$','Interpreter','latex')
legend("Optimal allocation w/o p# ", "Symmetric allocation w/ p#", "Optimal allocation w/ p#")
title("Heterogeneous case: N = 3, K = 2")