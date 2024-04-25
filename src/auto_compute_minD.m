function [download] = auto_compute_minD(N,K,leak_constraint)
% Given leakage constriant L1, L2, ... LN
% minimize D
% varaibles : p 
% D must satisfy 1<= D <= D* (D*: no leakage and 1: all leaked)

% calculate D*
D_star = 0;
for i = 1:K
    D_star = D_star + (1/N)^(i-1);
end

% get query table
[leakage_table, normal_table] = auto_create_PIRtable(N,K);
sz_leakage_table = size(leakage_table);

cvx_begin
    L = cvx(zeros(N,1));
    variables p(sz_leakage_table(1),K) D
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

    minimize (D)
    subject to 
        p >= 0;
        sum(p,1) == 1;
        N-sum(p(1:N+factorial(N),:),1) <= (N-1)*D;
        L <= leak_constraint;
cvx_end
download = D;

