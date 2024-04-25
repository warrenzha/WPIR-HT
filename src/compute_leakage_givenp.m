% calculate leakage given n, k, r and p

% assign N K r p
N = 3;
K = 2;
r = [0.1 0.3 0.6];
c = 0.02;
c1 = 0.01;
p = [0.1000 ,0.1000;
    0.0000 ,0.0000;
    0.0000 ,0.0000;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05;
    0.05 ,0.05]


% compute D 
D = (N-sum(p(1:N+factorial(N),1),1))/(N-1)

% find query table
[leakage_table, normal_table] = auto_create_PIRtable(N,K);
sz_leakage_table = size(leakage_table);

% compute leakage
L = zeros(N,1);

% compute leakage Li 
for i = 1:N
    search_table = leakage_table(:,i,:);
    count_map = zeros(sz_leakage_table(1),sz_leakage_table(3));
    for k = 1:K
        for t = 1:sz_leakage_table(1)
            now_searching = search_table(t,1,k);
            idx = find(ismember(search_table,now_searching));
            if count_map(idx(1)) == 0
                M = zeros(size(p));
                M(idx) = 1;
                a = max(sum(M.*p,1));
                L(i) = L(i)+a;
                count_map = count_map+M;
            end 
        end
    end
end
leakage = r*L