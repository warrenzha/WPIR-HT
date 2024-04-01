%% general N K case
N = 3;
K = 3;
F = K;
epsilon = 1;
dp = exp(epsilon);

% % get random key and query table
% [Rand,~,normal_query_table] = auto_create_PIRtable_RK(N,K);
% Rand_split = split(Rand,'');
% Rand_split = reshape(Rand_split(Rand_split ~= ""),F,[]);

f = [-1, zeros(1,F-1)];
A = [];
temp = 1;
for t = 1:F-1
    A(temp, t) = 1;
    A(temp, t+1) = -dp;
    temp = temp+1;
end
% for t = 2:F
%     A(temp, t-1) = -dp;
%     A(temp, t) = 1;
%     temp = temp+1;
% end
A

sjnum = zeros(size(f));
for j = 0:K-1
    sjnum(j+1) = nchoosek(K-1,j)*((N-1)^j);
end
sjnum
szA = size(A);
b = zeros(1,szA(1));
Aeq = factorial(N)*sjnum.*ones(1,F);
beq = 1;
lb = zeros(1,F);

% optimal p
[p,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb);

% 
p
fval
D = (1/(N-1))*(N-factorial(N)*p(1))
lambda.ineqlin
% lambda.eqlin
% lambda.lower