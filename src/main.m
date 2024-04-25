%% see what p is when min L 
r = [0.1 0.3 0.6];
[leak] = auto_compute_maxL(3,2,r,1.3);

%% see what p is when min D
r = [0.1 0.3 0.6];
L_constraints = [1.07,1.01,1]';
r*L_constraints
[download] = auto_compute_minD(3,2,L_constraints)

%% plot (L,D) graph of N=3 vs. N=2
clc
clear all 
close all 

r3 = [0.2 0.2 0.6];
r2 = [0.5 0.5];
D3 = 1:0.1:1.33;
D2 = 1:0.1:1.5;

leak3 = zeros(length(D3),1);
leak2 = zeros(length(D2),1);

for i = 1:length(D3)
    [leak3(i)] = auto_compute_maxL(3,2,r3,D3(i));
end
for i = 1:length(D2)
    [leak2(i)] = auto_compute_maxL(2,2,r2,D2(i));
end

plot(leak3,D3)
hold on 
plot(leak2,D2)
