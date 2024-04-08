function [download] = auto_compute_minD_DP_LPreduced(N,K,epsilon)
% Given leakage constriant L1, L2, ... LN
% minimize D
% varaibles : p 

% calculate D* (bound)
D_star = 0;
for i = 1:K
    D_star = D_star + (N*exp(epsilon))^(1-i);
end

D_upper = 1 + (N^(K-1)-1) / ((N-1)*(exp(epsilon)+N^(K-1)-1));

% Leakage setting
F = K;
dp = exp(epsilon);


f = [-1, zeros(1,F-1)];
A = [];
temp = 1;

for t = 1:F-1
    A(temp, t) = 1;
    A(temp, t+1) = -dp;
    temp = temp+1;
end

for t = 2:F
    A(temp, t-1) = -dp;
    A(temp, t) = 1;
    temp = temp+1;
end

sjnum = zeros(size(f));
for j = 0:K-1
    sjnum(j+1) = nchoosek(K-1,j)*((N-1)^j);
end

szA = size(A);
b = zeros(1,szA(1));
Aeq = factorial(N)*sjnum.*ones(1,F);
beq = 1;
lb = zeros(1,F);

% optimal p
[p,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb);

download = (1/(N-1))*(N-factorial(N)*p(1));