%% general N K case
N = 3;
K = 2;
F = N^(K-1);
epsilon = 1;
dp = exp(epsilon);

% get random key and query table
[Rand,~,normal_query_table] = auto_create_PIRtable_RK(N,K);
Rand_split = split(Rand,'');
Rand_split = reshape(Rand_split(Rand_split ~= ""),F,[]);

f = [factorial(N), (N/(N-1))*factorial(N)*ones(1,F-1)];
A = [];
temp = 1;
for t = 1:F
    count = zeros(F,1);
    for k = 1:K-1
        Rand_split_new = Rand_split;
        for i = 1:N
            Rand_split_new(t,k) = num2str(i-1);
            search_now = strjoin(Rand_split_new(t,:),'');
            count = count + ismember(Rand,search_now);
        end
    end
    indx = unique(find(count));
    for j = 1:length(indx)
        if t~= indx(j)
            A(temp, indx(j)) = -dp;
            A(temp, t) = 1;
            temp = temp+1;
        end
    end
end

szA = size(A);
b = zeros(1,szA(1));
Aeq = factorial(N)*ones(1,F);
beq = 1;
lb = zeros(1,F);

% optimal p
[p,fval,exitflag,output,lambda] = linprog(f,A,b,Aeq,beq,lb);

% 
p
fval
% lambda.ineqlin
% lambda.eqlin
% lambda.lower