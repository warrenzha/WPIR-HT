function [total_leakage,p] = auto_compute_MI_adjustp(N,K,r,D,mode,sharp)
% Given leakage weight r1 ... rN and download constriant D
% minimize MI
% varaibles : p & D
% D must satisfy 1<= D <= D* (D*: no leakage and 1: all leaked)
% sharp = 0: without p#

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
        % pw not uniform
        p(N+factorial(N)+1:end,:) = p1(N+2:end,:); 

    elseif mode == 2 % pf uniform (here w means every different w not cardinality)
        variables p1(N+factorial(N)+(N^(K-1)-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#, p0
        p(1:N+factorial(N),:) = p1(1:N+factorial(N),:);
        % pw uniform
        for f = 1:N^(K-1)-1
            for k = 1:K
                p(N+factorial(N)+factorial(N)*(f-1)+1:N+factorial(N)+factorial(N)*f,k) = p1(N+factorial(N)+f,k);
            end
        end
    
    elseif mode == 3 % paper allocation
        variables p1(N+N^(K-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#
        p(1:N,:) = p1(1:N,:);
        % p0 pf uniform
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
        % p0 pf uniform
        for f = 1:N^(K-1)
            for k = 1:K
                p(N+factorial(N)*(f-1)+1:N+factorial(N)*f,k) = p1(N+f,k);
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

    elseif mode == 6 % p# and p0 and p only
        variables p1(N+2,K)
        p = cvx(zeros(sz_leakage_table(1),K));
        %p#
        p(1:N,:) = p1(1:N,:);
        for k = 1:K
            p(N+1:N+factorial(N),k) = p1(N+1,k);
            p(N+factorial(N)+1:end,k) = p1(N+2,k);
        end
        
    elseif mode == 7 % p# and p0 and pw(w = cardinality f) only
        variables p1(N+K,K)
        p = cvx(zeros(sz_leakage_table(1),K));
        % p#
        p(1:N,:) = p1(1:N,:);
        % p0
        for k = 1:K
            p(N+1:N+factorial(N),k) = p1(N+1,k);
        end
        % pw
        wcount = N+factorial(N)+1;
        for w = 1:K-1
            for k = 1:K
                p(wcount:wcount+factorial(N)*(N-1)*N^(w-1)-1,k) = p1(N+w+1,k);
            end
            wcount = wcount+factorial(N)*(N-1)*N^(w-1);
        end

    elseif mode == 8 % p0 pw uniform
        variables p1(1+N^(K-1),K)
        p = cvx(zeros(sz_leakage_table(1),K));

        %p#
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
    
    I = cvx(zeros(N,1));
    % compute leakage Li 
    for i = 1:N
        search_table = leakage_table(:,i,:);
        for m = 1:K
            count_map = zeros(sz_leakage_table(1),sz_leakage_table(3));
            %for k = 1:K
                for t = 1:sz_leakage_table(1)
                    now_searching = search_table(t,1,m);
                    idx = find(ismember(search_table,now_searching));
                    if count_map(idx(1)) == 0
                        M = zeros(size(p));
                        M(idx) = 1;
                        a = sum(M.*p,1);
                        x = a(m);
                        y = sum(a)*(1/K);
                        I(i) = I(i)+ (1/K)*rel_entr(x,y);
                        count_map = count_map+M;
                    end 
                end
            %end
        end
    end

    minimize (r*I)
    subject to 
        if sharp == 0
            p(1:N,:) == 0;
        end
        p >= 0;
        sum(p,1) == 1;
        N-sum(p(1:N+factorial(N),:),1) <= (N-1)*D;
cvx_end
total_leakage = r*I;
I
p