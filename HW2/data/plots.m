clc;clear;close all;

D = csvread("../T.csv");

xs = D(1,:);
ys = D(2,:);
D(1:2,:) = [];
N = length(xs);


[X, Y] =  meshgrid(xs, ys);

Q = zeros(N, N);
S = zeros(N, N);

for i = 1:N
    for j = 1:N
        Q(i,j) = 4*ys(j)^3 - 6*ys(j)^2 + 2 - 6*(1 - xs(i)^2)*(2*ys(j) - 1);
        S(i,j) = (1 - xs(i)^2)*(2*ys(j)^3 - 3*ys(j)^2 + 1);
    end
end

surf(X, Y, D'-S',"edgecolor","blue","facecolor","blue","facealpha",1);hold on;
% surf(X, Y, Q')
% surf(X, Y, S',"edgecolor","green","facecolor","green","facealpha",1);
xlabel("x")
ylabel("y")
