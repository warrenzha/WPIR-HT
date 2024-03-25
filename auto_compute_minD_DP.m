function [download] = auto_compute_minD_DP(N,K,epsilon)
% Given leakage constriant L1, L2, ... LN
% minimize D
% varaibles : p 
% D must satisfy 1<= D <= D* (D*: no leakage and 1: all leaked)

% calculate D*
D_star = 0;
for i = 1:K
    D_star = D_star + (N*exp(epsilon))^(1-i);
end


% get query table
[leakage_table, normal_table] = auto_create_PIRtable(N,K);
sz_normal_table = size(normal_table);
q_numb_normal = N^K;

cvx_begin
    L = cvx(zeros(N,1));
    same_q_prob = cvx(zeros(q_numb_normal,K,N));
    variables p(sz_normal_table(1),K) D
    % compute leakage Li 
    for i = 1:N
        search_table = normal_table(:,i,:);
        count_map = zeros(sz_normal_table(1),sz_normal_table(3));
        count_q = 0;
        for k = 1:K
            for t = 1:sz_normal_table(1)
                now_searching = search_table(t,1,k);
                idx = find(ismember(search_table,now_searching));
                if count_map(idx(1)) == 0
                    M = zeros(size(p));
                    M(idx) = 1;
                    a = max(sum(M.*p,1));
                    L(i) = L(i)+a;
                    count_map = count_map+M;

                    count_q = count_q + 1;
                    same_q_prob(count_q,:,i) = sum(M.*p,1);
                end 
            end
        end
    end

    minimize (D)
    subject to 
        p >= 0;
        sum(p,1) == 1;
        N-sum(p(1:factorial(N),:),1) == (N-1)*D;
        for k1 = 1:K
            for k2 = 1:K
                if k1 ~= k2
                    same_q_prob(:,k1,:) <= exp(epsilon).*same_q_prob(:,k2,:);
                end
            end
        end
cvx_end
download = D;
p
D_star
