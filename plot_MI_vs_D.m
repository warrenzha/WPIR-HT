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

% D = 1.00:0.02:D_star;
D = linspace(1.00,D_star,17);

% create initial matirx 
MI_leakage_homo_1 = zeros(size(D));
MI_leakage_homo_2 = zeros(size(D));
MI_leakage_homo_3 = zeros(size(D));
MI_leakage_homo_4 = zeros(size(D));

MI_leakage_hete_1 = zeros(size(D));
MI_leakage_hete_2 = zeros(size(D));
MI_leakage_hete_3 = zeros(size(D));
MI_leakage_hete_4 = zeros(size(D));

%% Homo
for d = 1:length(D)
    % line 1 : weight same optimum without #
    [MI_leakage_homo_1(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_same,D(d),0,0);
    
    % line 2 : weight same symmetric allocation with #
    [MI_leakage_homo_2(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_same,D(d),5,1);

    % line 3 : weight same paper allocation we found with # 
    [MI_leakage_homo_3(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_same,D(d),4,1);
    
    % line 3 : weight same paper allocation we found with # 
    [MI_leakage_homo_4(d), p_homo_rsame_opt] = auto_compute_MI_adjustp(N,K,r_same,D(d),0,1);
end

%% Hetero
for d = 1:length(D)
    % line 1 : weight diff optimum without #
    [MI_leakage_hete_1(d), p_hete_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),0,0);

    % line 2 : weight diff symmetric allocation with #
    [MI_leakage_hete_2(d), p_hete_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),5,1);
    
    % line 3 : weight diff paper allocation we found with # 
    [MI_leakage_hete_3(d), p_hete_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),4,1);

    % line 4 : weight diff optimum with # 
    [MI_leakage_hete_4(d), p_hete_rsame_opt] = auto_compute_MI_adjustp(N,K,r_diff,D(d),4,1);
end

%% plot homo
figure
set(gcf,'position',[50, 100, 550, 400])
set(gcf,'Color',[0.9,0.9,0.9])
plot(MI_leakage_homo_1,D,'-square','LineWidth',1.2)
hold on
plot(MI_leakage_homo_2,D,'-o','LineWidth',1.2)
hold on
plot(MI_leakage_homo_3,D,'-*','LineWidth',1.2)
hold on
plot(MI_leakage_homo_4,D,'-.','LineWidth',1.2)

grid on
xlim([0 0.3])
xlabel('MI constraint $\rho$','Interpreter','latex')
ylabel('Download cost $D$','Interpreter','latex')
legend("TSC (w/o p#)", "Symmetric WPIR#", "Optimal WPIR#", "Numerical solution")
title("Homogeneous case: N = 3, K = 2")

%% plot hetero
figure
set(gcf,'position',[50, 100, 550, 400])
set(gcf,'Color',[0.9,0.9,0.9])
plot(MI_leakage_hete_1,D,'-square','LineWidth',1.2)
hold on
plot(MI_leakage_hete_2,D,'-o','LineWidth',1.2)
hold on
plot(MI_leakage_hete_3,D,'-*','LineWidth',1.2)
hold on
plot(MI_leakage_hete_4,D,'-.','LineWidth',1.2)

grid on
xlim([0 0.08])
xlabel('MI constraint $\rho$','Interpreter','latex')
ylabel('Download cost $D$','Interpreter','latex')
legend("TSC (w/o p#)", "Symmetric WPIR#", "Optimal WPIR#", "Numerical solution")
title("Heterogeneous case: N = 3, K = 2")