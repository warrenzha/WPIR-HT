function [leakage_query_table,normal_query_table] = auto_create_PIRtable(N,K)
% servers N and messages K 
% do permutation of 0:N-1 
Per = perms(0:N-1);

% set random key F
cardi_F = N^(K-1);
F = strings(cardi_F,1);
for i = 1:cardi_F
    if strlength(dec2base(i-1,N)) < (K-1)
        gap = (K-1)-strlength(dec2base(i-1,N));
        F(i,1) = strjoin(repmat("0", 1, gap),'')+dec2base(i-1,N);
    else
        F(i,1) = dec2base(i-1,N);
    end
end

% create PIR table 
normal_query_table = strings(length(Per)*length(F), N, K);
for k = 1:K
    for n = 0:N-1
        i = 0;
        for f = 1:length(F)
            for per = 1:length(Per)
                i = i+1;
                % sum random key f
                total_randomf = split(F(f,1),"");
                total_randomf = sum(str2double(total_randomf(strlength(total_randomf)>0)));

                % compute query
                per_n = mod(Per(per,n+1)-total_randomf,N);
                q = insertAfter(F(f,1),k-1,string(per_n));
                normal_query_table(i,n+1,k) = q;
            end
        end
    end
end

% create leakage table 
leakage_query_table = strings(length(Per)*length(F)+N, N, K);
for i = 1:N
    for n = 1:N
        for k = 1:K
            if i == n
                leakage_query_table(i,n,k) = "#"+string(k);
            else 
                leakage_query_table(i,n,k) = strjoin(repmat("0", 1, K),'');
            end
        end
    end
end
leakage_query_table(N+1:end,:,:) = normal_query_table;
