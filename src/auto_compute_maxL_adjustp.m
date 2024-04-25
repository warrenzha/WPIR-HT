function [total_leakage,p] = auto_compute_maxL_adjustp(N,K,r,D, mode,sharp)
% Given leakage weight r1 ... rN and download constriant D
% minimize MaxL 
% varaibles : p & D
% D must satisfy 1<= D <= D* (D*: no leakage and 1: all leaked)
% sharp = 0: without p#
% mode = 1: p0 uniform/ mode = 2: pw uniform / mode = 3: p0,pw uniform
% mode = 4: paper allocation/ mode 5 : symmetric allocation

% calculate D*
D_star = 0;
for i = 1:K
    D_star = D_star + (1/N)^(i-1);
end

% get query table
[leakage_table, normal_table] = auto_create_PIRtable(N,K);
sz_leakage_table = size(leakage_table);

% begin convex problem
cvx_begin
    % mode
    if mode == 1 %p0 uniform
        variables p1(N+1+factorial(N)*(N^(K-1)-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#
        p(1:N,:) = p1(1:N,:);
        % p0 uniform part
        for k = 1:K 
            p(N+1:N+factorial(N),k) = p1(N+1,k);
        end
        % pf not uniform
        p(N+factorial(N)+1:end,:) = p1(N+2:end,:); 

    elseif mode == 2 % pf uniform
        variables p1(N+factorial(N)+(N^(K-1)-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#, p0
        p(1:N+factorial(N),:) = p1(1:N+factorial(N),:);
        % pf uniform
        for f = 1:N^(K-1)-1
            for k = 1:K
                p(N+factorial(N)+factorial(N)*(f-1)+1:N+factorial(N)+factorial(N)*f,k) = p1(N+factorial(N)+f,k);
            end
        end
    
    elseif mode == 3 % p0 pf uniform
        variables p1(N+N^(K-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#
        p(1:N,:) = p1(1:N,:);
        % p0 pw uniform
        for f = 1:N^(K-1)
            for k = 1:K
                p(N+factorial(N)*(f-1)+1:N+factorial(N)*f,k) = p1(N+f,k);
            end
        end
        
    elseif mode == 4 % paper allocation
        variables p1(N+N^(K-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#
        p(1:N,:) = p1(1:N,:);
        % p0 pw
        for f = 1:N^(K-1)
            for k = 1:K
                for n = 1:N
                    p(N+factorial(N)*(f-1)+1+(n-1)*factorial(N-1),k) = p1(N+f,k);
                end
            end
        end
        
    elseif mode == 5 % p0 pw uniform, p# symmetric
        variables p1(1+N^(K-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#
        for f = 1:N
            p(f,:) = p1(1,:);
        end
        % p0 pw uniform
        for f = 1:N^(K-1)
            for k = 1:K
                p(N+factorial(N)*(f-1)+1:N+factorial(N)*f,k) = p1(1+f,k);
            end
        end
        
    else 
        variables p(sz_leakage_table(1),K)
    end
    
    % compute leakage L(1:N)
    L = cvx(zeros(N,1));
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

    minimize (r*L)
    subject to 
        if sharp == 0
            p(1:N,:) == 0;
        end
        p >= 0;
        sum(p,1) == 1;
        N-sum(p(1:N+factorial(N),:),1) <= (N-1)*D;
cvx_end

total_leakage = r*L;
L
p