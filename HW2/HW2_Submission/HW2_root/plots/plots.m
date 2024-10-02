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

surf(X, Y, D',"edgecolor","black","facecolor","blue","facealpha",1);hold on;
% surf(X, Y, Q')
% surf(X, Y, S',"edgecolor","green","facecolor","green","facealpha",1);
xlabel("x")
ylabel("y")

%%
figure;
res = csvread("residuals_c_j_10.csv");
semilogy(res);

%% Part c
clc;clear;close all;

D10 = csvread("T_c_J_10.csv");
D20 = csvread("T_c_J_20.csv");
D40 = csvread("T_c_J_40.csv");

xs10 = D10(1,:);
ys10 = D10(2,:);
D10(1:2,:) = [];
N10 = length(xs10);

xs20 = D20(1,:);
ys20 = D20(2,:);
D20(1:2,:) = [];
N20 = length(xs20);

xs40 = D40(1,:);
ys40 = D40(2,:);
D40(1:2,:) = [];
N40 = length(xs40);

S = zeros(N40, N40);

for i = 1:N40
    for j = 1:N40
        S(i,j) = (1 - xs40(i)^2)*(2*ys40(j)^3 - 3*ys40(j)^2 + 1);
    end
end

[X10, Y10] = meshgrid(xs10, ys10);
[X20, Y20] = meshgrid(xs20, ys20);
[X40, Y40] = meshgrid(xs40, ys40);

error10 = sum(sum((D10 - S(1:4:end, 1:4:end)).^2))./11^2;
error20 = sum(sum((D20 - S(1:2:end, 1:2:end)).^2))./21^2;
error40 = sum(sum((D40 - S(:,:)).^2))./41^2;

figure("position",[50 50 1500 400])
subplot(1,3,1)
contour(X10, Y10, D10,'r');hold on;
contour(X20, Y20, D20,'g');
contour(X40, Y40, D40,'b');
contour(X40, Y40, S,'k');
axis equal
xlabel("X Coordinate")
ylabel("Y Coordinate")
title("Contours at Different Resolution")
legend(["N = 10","N = 20","N = 40","Exact"],'location','northwest')
fontsize(gca,16,"points")

subplot(1,3,2)
plot(ys10, D10(6,:),'r'); hold on;
plot(ys20, D20(11,:),'g')
plot(ys40, D40(21,:),'b')
plot(ys40, S(21,:),'k')
xlabel("Y Coordinate");
ylabel("Temperature");
title("Temperature along T(0.5,y)");
legend(["N = 10","N = 20","N = 40","Exact"],'location','southwest')
fontsize(gca, 16, "points")

N = [10 20 40];
N2 = [12, 35];

subplot(1,3,3)
loglog(N, [error10, error20, error40],'k-','linewidth',2);hold on;
loglog(N, [error10, error20, error40],'k.','markersize',20);
loglog(N2, 1e-1*(1./N2).^4,'k-')
text(20,0.1e-5,"slope = -4","rotation",-40);
title("L_2 error");
xlabel("N")
ylabel("L_2 error");
fontsize(gca,16,"points");

%% Part D

clc;clear;close all;

D10 = csvread("T_d_GS_10.csv");
D20 = csvread("T_d_GS_20.csv");
D40 = csvread("T_d_GS_40.csv");
D80 = csvread("T_d_GS_80.csv");
D160 = csvread("T_d_GS_160.csv");

xs10 = D10(1,:);
ys10 = D10(2,:);
D10(1:2,:) = [];
N10 = length(xs10);

xs20 = D20(1,:);
ys20 = D20(2,:);
D20(1:2,:) = [];
N20 = length(xs20);

xs40 = D40(1,:);
ys40 = D40(2,:);
D40(1:2,:) = [];
N40 = length(xs40);

xs80 = D80(1,:);
ys80 = D80(2,:);
D80(1:2,:) = [];
N80 = length(xs80);

xs160 = D160(1,:);
ys160 = D160(2,:);
D160(1:2,:) = [];
N160 = length(xs160);

S = zeros(N160, N160);

for i = 1:N160
    for j = 1:N160
        S(i,j) = (1 - xs160(i)^2)*(2*ys160(j)^3 - 3*ys160(j)^2 + 1);
    end
end

[X10, Y10] = meshgrid(xs10, ys10);
[X20, Y20] = meshgrid(xs20, ys20);
[X40, Y40] = meshgrid(xs40, ys40);
[X80, Y80] = meshgrid(xs80, ys80);
[X160, Y160] = meshgrid(xs160, ys160);

error10 = sum(sum((D10 - S(1:16:end, 1:16:end)).^2))./11^2;
error20 = sum(sum((D20 - S(1:8:end, 1:8:end)).^2))./21^2;
error40 = sum(sum((D40 - S(1:4:end, 1:4:end)).^2))./41^2;
error80 = sum(sum((D80 - S(1:2:end, 1:2:end)).^2))./81^2;
error160 = sum(sum((D160 - S(:,:)).^2))./161^2;

N = [10 20 40 80 160];
E = [error10, error20, error40, error80, error160];
N2 = [15, 100];

loglog(N, E, "k-",'linewidth',2);hold on;
loglog(N, E, "k.",'markersize',25)
xlabel("N cells");
ylabel("L_2 error")
title("L_2 error convergence for Gauss-Seidel");
text(40,6e-8,"slope = -4","rotation",-43);
loglog(N2, 1e-1*(1./N2).^4,'k-')
ylim([5e-10,1e-5])
fontsize(gca,16,"points");
 
%% Part E

clc;clear;close all;

DJ = csvread("residuals_e_J.csv");
DGS = csvread("residuals_e_GS.csv");
DSOR11 = csvread("residuals_e_SOR11.csv");
DSOR12 = csvread("residuals_e_SOR12.csv");
DSOR13 = csvread("residuals_e_SOR13.csv");
DSOR14 = csvread("residuals_e_SOR14.csv");
DSOR15 = csvread("residuals_e_SOR15.csv");
DSOR16 = csvread("residuals_e_SOR16.csv");
DSOR17 = csvread("residuals_e_SOR17.csv");
DSOR18 = csvread("residuals_e_SOR18.csv");
DSOR19 = csvread("residuals_e_SOR19.csv");

NJ = 1:1:length(DJ);
NGS = 1:1:length(DGS);
NSOR11 = 1:1:length(DSOR11);
NSOR12 = 1:1:length(DSOR12);
NSOR13 = 1:1:length(DSOR13);
NSOR14 = 1:1:length(DSOR14);
NSOR15 = 1:1:length(DSOR15);
NSOR16 = 1:1:length(DSOR16);
NSOR17 = 1:1:length(DSOR17);
NSOR18 = 1:1:length(DSOR18);
NSOR19 = 1:1:length(DSOR19);

semilogy(NJ, DJ, 'r','linewidth',2);hold on;
semilogy(NGS, DGS, 'g','linewidth',2);hold on;
semilogy(NSOR11, DSOR11,'b-','linewidth',2);
semilogy(NSOR12, DSOR12,'b-.','linewidth',2);
semilogy(NSOR13, DSOR13,'b--','linewidth',2);
semilogy(NSOR14, DSOR14,'m-','linewidth',2);
semilogy(NSOR15, DSOR15,'m-.','linewidth',2);
semilogy(NSOR16, DSOR16,'m--','linewidth',2);
semilogy(NSOR17, DSOR17,'c-','linewidth',2);
semilogy(NSOR18, DSOR18,'c-.','linewidth',2);
semilogy(NSOR19, DSOR19,'c--','linewidth',2); 
legend(["Jacobi","Gauss-Seidel","SOR \omega=1.1", "SOR \omega=1.2","SOR \omega=1.3", "SOR \omega=1.4", ...
    "SOR \omega=1.5", "SOR \omega=1.6", "SOR \omega=1.7", "SOR \omega=1.8", ...
    "SOR \omega=1.9"],"location","eastoutside");
xlabel("Iteration Count");
ylabel("Residual");
title("Residual vs. Iteration");
fontsize(gca,16,"points");

%% Part f

clc;clear;close all;

N = [10 20 40 80];
N2 = [15 50];
tSolveJ = [1.4e-3 2.0e-3, 2.9e-2, 0.47];
tIterJ = [3.2e-7 1.2e-6, 4.2e-6, 1.7e-5];
tSolveGS = [1.85e-4, 1.8e-4, 4.5e-2, 0.73];
tIterGS = [8.5e-7, 8.4e-7, 1.3e-5, 5.25e-5];
tSolveSOR = [4.1e-4, 6.5e-4, 1.7e-3, 2.6e-2];
tIterSOR = [6.1e-7, 2.5e-6, 8.7e-6, 3.76e-5];

figure("position",[50 50 1500 400])
subplot(1,3,1)
loglog(N,tSolveJ,'b-','linewidth',2);hold on;
loglog(N, tIterJ,'r-','linewidth',2);
plot(N2, 1e-8*N2.^2,'k-')
text(25,1e-5,"slope = 2","rotation", 13);
plot(N2, 5e-8*N2.^4,'k-')
text(25, 4e-2,"slope = 4","rotation", 26);
xlabel("N Points");
ylabel("Runtime (s)");
title("Jacobi Runtime");
fontsize(gca,16,"points");

subplot(1,3,2)
loglog(N,tSolveGS,'b-','linewidth',2);hold on;
loglog(N, tIterGS,'r-','linewidth',2);
plot(N2, 1e-8*N2.^2,'k-')
text(25,1e-5,"slope = 2","rotation", 13);
plot(N2, 5e-8*N2.^4,'k-')
text(25, 4e-2,"slope = 4","rotation", 26);
ylim([1e-7, 1])
xlabel("N Points");
ylabel("Runtime (s)");
title("Gauss-Seidel Runtime");
fontsize(gca,16,"points");

subplot(1,3,3)
loglog(N,tSolveSOR,'b-','linewidth',2);hold on;
loglog(N, tIterSOR,'r-','linewidth',2);
plot(N2, 1e-8*N2.^2,'k-')
text(25,1e-5,"slope = 2","rotation", 15);
plot(N2, 2e-9*N2.^4,'k-')
text(25, 1.5e-3,"slope = 4","rotation", 26);
ylim([1e-7, 1])
xlabel("N Points");
ylabel("Runtime (s)");
title("SOR \omega=1.9 Runtime");
fontsize(gca,16,"points");

%% Part g
clc;clear;close all;

N = [10 20 40 80];
N2 = [15 50];

tSolveGS = [1.85e-4, 1.8e-4, 4.5e-2, 0.73];
tSolveMG = [2.65e-6, 9.7e-6, 3.6e-5, 1.83e-4];
tSolveJ = [1.4e-3 2.0e-3, 2.9e-2, 0.47];
tSolveSOR = [4.1e-4, 6.5e-4, 1.7e-3, 2.6e-2];

loglog(N, tSolveJ, 'g','linewidth',2);hold on;
loglog(N, tSolveGS,'r-','linewidth',2);
loglog(N, tSolveSOR,'m-','linewidth',2);
loglog(N, tSolveMG, 'b-','linewidth',2);
legend(["Jacobi", "Gauss-seidel", "SOR", "Algebraic Multigrid", "slope = 2", "slope = 4"], 'location','northwest')


plot(N2, 1e-7*N2.^2,'k-')
plot(N2, 4e-8*N2.^4,'k--')
xlabel("N");
ylabel("Time to convergence")
fontsize(gca,16,"points")