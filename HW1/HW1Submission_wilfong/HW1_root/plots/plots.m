%% Part b

clc;clear;close all;

D10 = readtable("outputB10.csv");
D20 = readtable("outputB20.csv");
D30 = readtable("outputB30.csv");
D40 = readtable("outputB40.csv");
D50 = readtable("outputB50.csv");

xs = linspace(0, 1, 1000);

U = 1;
G = 0.1;
Q = 0;
L = 1;

ref = @(x) 1 - (exp(x.*U/G) - 1)/(exp(L*U/G) - 1);

E10 = mean(D10.Var2 - ref(D10.Var1));
E20 = mean(D20.Var2 - ref(D20.Var1));
E30 = mean(D30.Var2 - ref(D30.Var1));
E40 = mean(D40.Var2 - ref(D40.Var1));
E50 = mean(D50.Var2 - ref(D50.Var1));

E = [E10 E20 E30 E40 E50];
N = [10 20 30 40 50];

T = [0.69 1.16 1.64 2.10 2.56];

figure('position',[50 50 1000 300]);
subplot(1,3,1);
plot(D10.Var1, D10.Var2,'linewidth',2);hold on;
plot(D20.Var1, D20.Var2,'linewidth',2);hold on;
plot(D30.Var1, D30.Var2,'linewidth',2);hold on;
plot(D40.Var1, D40.Var2,'linewidth',2);hold on;
plot(D50.Var1, D50.Var2,'linewidth',2);hold on;
plot(xs, ref(xs),'-k','linewidth',1)
legend("N = 10", "N - 20", "N = 30", "N = 40", "N = 50", "Ref.", "location","southwest");
xlabel("x");
ylabel("f(x)");
title("f(x) vs. x");
fontsize(gca,16,"points");

subplot(1,3,2);
loglog(N, E,'k--');hold on;
loglog(N, E, 'k.','markersize',15)
loglog(N, (1./N).^2,'k-')
xlabel("N");
ylabel("Average Error");
title("Average Error vs. N");
text(20, 3e-3, "Slope = -2","rotation",-42)
fontsize(gca, 16,"points");

subplot(1,3,3);
loglog(N, T,'k--');hold on;
loglog(N, T, 'k.','markersize',15)
loglog(N, 0.1.*N.^1,'k-')
xlabel("N");
ylabel("Runtime (\mu s)");
title("Runtime vs. N");
ylim([0.5 6])
yticks([0.5 1 2 3 4 5 6]);
text(18, 2.1, "Slope = 1","rotation", 38)
fontsize(gca, 16,"points");
%% Part C
clc;clear;close all;

D1 = readtable("outputC1.csv");
D2 = readtable("outputC2.csv");
D3 = readtable("outputC3.csv");

ref1 = @(x) -1 + sqrt(4 - 3.*x);
ref2 = @(x) -1 + sqrt(4 - 2.*x - x.^2);
ref3 = @(x) -1 + sqrt(4 - (8/3).*x - (x.^3)/3);

figure('position',[50 50 700 300]);
subplot(1,2,1);
plot(D1.Var1, D1.Var2,'linewidth',2);hold on;
plot(D2.Var1, D2.Var2,'linewidth',2);hold on;
plot(D3.Var1, D3.Var2,'linewidth',2);hold on; 
legend({'$\varphi = 0$', '$\varphi = 0.1$', '$\varphi = 0.1x$'},'Interpreter','latex', 'location','southwest');
xlabel("x");
ylabel("f(x)");
title("f(x) vs. x");
fontsize(gca,16,"points");

subplot(1,2,2);
plot(D1.Var1, abs(D1.Var2 - ref1(D1.Var1)),'linewidth',2);hold on;
plot(D2.Var1, abs(D2.Var2 - ref2(D2.Var1)),'linewidth',2);hold on;
plot(D3.Var1, abs(D3.Var2 - ref3(D3.Var1)),'linewidth',2);hold on;
xlabel("x");
ylabel("|error|");
title("|error| vs. x");
fontsize(gca, 16,"points");

