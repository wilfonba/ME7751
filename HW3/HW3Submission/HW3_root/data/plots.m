clc;close all;clear;

D = csvread("../output.csv");

figure("Position",[50 50 1600 300])
for i = 2:size(D,1)
    cla;
    plot(D(1,:), D(i,:),'k-','linewidth',2);
    ylim([-1 2])
    drawnow;
end

%% Problem 2
clc;clear;close all;

C1Ref = @(x, t) ((4.*pi.*0.01.*t).^(-0.5)).*(exp((-(x - t).^2)/(4.*0.01.*t)));
C2Ref = @(x, t) (0.4*pi).^(-0.5).*exp(-2.5.*(x - t).^2);

dt = 0.1;

figure("position",[50 50 1000 1000]);
% Case 1
D0 = csvread("dataQ2C1TS0.csv");
D1 = csvread("dataQ2C1TS1.csv");
D2 = csvread("dataQ2C1TS2.csv");


subplot(4,5,1)
plot(D0(1,:), D0(3,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(3,:),'g-','linewidth',1)
plot(D2(1,:), D2(3,:),'b-','linewidth',1)
plot(D0(1,:), C1Ref(D0(1,:),20),'k-','linewidth',1)
xlim([19, 21])
ylim([0.3, 0.7])
set(gca, 'XTick', []);
set(gca, 'YTick', []);

subplot(4,5,2:3)
plot(D0(1,:), D0(3,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(3,:),'g-','linewidth',1)
plot(D2(1,:), D2(3,:),'b-','linewidth',1)
plot(D0(1,:), C1Ref(D0(1,:),20),'k-','linewidth',1)
plot([19 19 21 21 19],[0.3 0.7 0.7 0.3 0.3],'k-');
xlim([5 45])
ylim([-0.1 0.8])
xlabel("x-position")
ylabel("\phi")
title("Case 1, T = 20")

subplot(4,5,6)
plot(D0(1,:), D0(4,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(4,:),'g-','linewidth',1)
plot(D2(1,:), D2(4,:),'b-','linewidth',1)
plot(D0(1,:), C1Ref(D0(1,:),30),'k-','linewidth',1)
xlim([28.5, 31.5])
ylim([0.2, 0.6])
set(gca, 'XTick', []);
set(gca, 'YTick', []);

subplot(4,5,7:8)
plot(D0(1,:), D0(4,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(4,:),'g-','linewidth',1)
plot(D2(1,:), D2(4,:),'b-','linewidth',1)
plot(D0(1,:), C1Ref(D0(1,:),30),'k-','linewidth',1)
plot([28.5 28.5 31.5 31.5 28.5],[0.2 0.6 0.6 0.2 0.2],'k-')
xlim([5 45])
ylim([-0.1 0.8])
xlabel("x-position")
ylabel("\phi")
title("Case 1, T = 30")

subplot(4,5,11)
plot(D0(1,:), D0(5,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(5,:),'g-','linewidth',1)
plot(D2(1,:), D2(5,:),'b-','linewidth',1)
plot(D0(1,:), C1Ref(D0(1,:),40),'k-','linewidth',1)
xlim([38, 42])
ylim([0.1, 0.5])
set(gca, 'XTick', []);
set(gca, 'YTick', []);

subplot(4,5,12:13)
plot(D0(1,:), D0(5,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(5,:),'g-','linewidth',1)
plot(D2(1,:), D2(5,:),'b-','linewidth',1)
plot(D0(1,:), C1Ref(D0(1,:),40),'k-','linewidth',1)
plot([38 38 42 42 38],[0.1 0.5 0.5 0.1 0.1],'k-')
xlim([5 45])
ylim([-0.1 0.8])
xlabel("x-position")
ylabel("\phi")
title("Case 1, T = 40")

subplot(4, 5, 16:20)
plot(D0(1,:) + 100, D0(5,:),'r-','linewidth',2);hold on;
plot(D1(1,:) + 100, D1(3,:),'g-','linewidth',2)
plot(D2(1,:) + 100, D2(3,:),'b-','linewidth',2)
plot(D2(1,:) + 100, D2(3,:),'k-','linewidth',2)
plot(D2(1,:), D2(2,:),'b-');
plot(D2(1,:),D2(3,:),'b-');
plot(D2(1,:),D2(4,:),'b-');
plot(D2(1,:),D2(5,:),'b-');
text(9, 1,"T = 10")
text(19, 0.76,"T = 20")
text(29, 0.65,"T = 30")
text(39, 0.59,"T = 40")
xlabel("x-positiona")
ylabel("\phi")
xlim([5 45])
ylim([-0.1, 1.1])
title("Case 1 Crank-Nicolson")
legend(["Explicit Euler","Implicit Euler","Crank-Nicolson","Reference"],'location','westoutside','fontsize',16)

% Case 2
D0 = csvread("dataQ2C2TS0.csv");
D1 = csvread("dataQ2C2TS1.csv");
D2 = csvread("dataQ2C2TS2.csv");

subplot(4,5,4:5)
plot(D0(1,:), D0(3,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(3,:),'g-','linewidth',1)
plot(D2(1,:), D2(3,:),'b-','linewidth',1)
plot(D0(1,:), C2Ref(D0(1,:),20),'k-','linewidth',1)
xlim([5 45])
ylim([-0.3 1])
xlabel("x-position")
ylabel("\phi")
title("Case 2, T = 20")

subplot(4,5,9:10)
plot(D0(1,:), D0(4,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(4,:),'g-','linewidth',1)
plot(D2(1,:), D2(4,:),'b-','linewidth',1)
plot(D0(1,:), C2Ref(D0(1,:),30),'k-','linewidth',1)
xlim([5 45])
ylim([-0.3 1])
xlabel("x-position")
ylabel("\phi")
title("Case 2, T = 30")

subplot(4,5,14:15)
plot(D0(1,:), D0(5,:),'r-','linewidth',1);hold on;
plot(D1(1,:), D1(5,:),'g-','linewidth',1)
plot(D2(1,:), D2(5,:),'b-','linewidth',1)
plot(D0(1,:), C2Ref(D0(1,:),40),'k-','linewidth',1)
xlim([5 45])
ylim([-0.3 1])
xlabel("x-position")
ylabel("\phi")
title("Case 2, T = 40")

%% Problem 3
clc;clear all;close all

D0 = csvread("dataQ3C1TS0.csv");
D1 = csvread("dataQ3C1TS1.csv");
D2 = csvread("dataQ3C1TS2.csv");

T = 10:0.1:40;

C1Ref = @(x, t) ((4.*pi.*0.01.*t).^(-0.5)).*(exp((-(x - t).^2)/(4.*0.01.*t)));

RefSol = zeros(length(T));

for i = 1:length(T)
    RefSol(i,1) = C1Ref(15, T(i));
end 

figure("position",[50 50 1000 250]) 
subplot(1, 2, 1);
plot(T, D0(2:end,126),'r-','linewidth',2);hold on;
plot(T, D1(2:end,126),'g-','linewidth',2);
plot(T, D2(2:end,126),'b-','linewidth',2);
plot(T, RefSol(:,1),'k-','linewidth',2);
legend(["Explicit Euler","Implicit Euler","Crank-Nicolson","Reference"])
title("x = 15");
ylabel("\phi(t)");
xlabel("t")
ylim([-0.1 0.9])
fontsize(gca, 16,"points")

for i = 1:length(T)
    RefSol(i,1) = C1Ref(25, T(i));
end 

subplot(1, 2, 2);
plot(T, D0(2:end,251),'r-','linewidth',2);hold on;
plot(T, D1(2:end,251),'g-','linewidth',2);
plot(T, D2(2:end,251),'b-','linewidth',2);
plot(T, RefSol(:,1),'k-','linewidth',2);
ylim([-0.1 0.9])
title("x = 25");
ylabel("\phi(t)");
xlabel("t")
fontsize(gca, 16,"points")

%% Problem 4

clc;clear all;close all

D0 = csvread("dataQ4N256.csv");
D1 = csvread("dataQ4N512.csv");
D2 = csvread("dataQ4N1024.csv");
D3 = csvread("dataQ4N2048.csv");
D4 = csvread("dataQ4N4096.csv");

N = [256, 512, 1024, 2048, 4096];
E = zeros(5,3);

C1Ref = @(x, t) ((4.*pi.*0.01.*t).^(-0.5)).*(exp((-(x - t).^2)/(4.*0.01.*t)));

E(1) = sqrt(sum((D0(end,:) - C1Ref(D0(1,:),40)).^2)/256);
E(2) = sqrt(sum((D1(end,:) - C1Ref(D1(1,:),40)).^2)/512);
E(3) = sqrt(sum((D2(end,:) - C1Ref(D2(1,:),40)).^2)/1024);
E(4) = sqrt(sum((D3(end,:) - C1Ref(D3(1,:),40)).^2)/2048);
E(5) = sqrt(sum((D4(end,:) - C1Ref(D4(1,:),40)).^2)/4096);

loglog(N, E,'k-','linewidth',2);hold on;
loglog(N, E,'k.','markersize',20);
loglog(N(2:4),2.5e3./(N(2:4).^2),'k-');
text(2^10, 3.5e-3,"Slope = -2","rotation",-45)
xlabel("Number of Grid Points")
xticks(N)
xlim([2^7, 2^13])
ylim([1e-4 3e-2])
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('2^{%d}', log2(x)), xticks, 'UniformOutput', false));
ylabel("L_2 error");
title("Spatial Convergence with Crank-Nicolson");
fontsize(gca,16,"points");

%% Problem 5
clc;clear all;close all;

D01 = csvread("dataQ5T1TS0.csv");
D02 = csvread("dataQ5T2TS0.csv");
D03 = csvread("dataQ5T3TS0.csv");
D04 = csvread("dataQ5T4TS0.csv");
D05 = csvread("dataQ5T5TS0.csv");

D11 = csvread("dataQ5T1TS1.csv");
D12 = csvread("dataQ5T2TS1.csv");
D13 = csvread("dataQ5T3TS1.csv");
D14 = csvread("dataQ5T4TS1.csv");
D15 = csvread("dataQ5T5TS1.csv");

D21 = csvread("dataQ5T1TS2.csv");
D22 = csvread("dataQ5T2TS2.csv");
D23 = csvread("dataQ5T3TS2.csv");
D24 = csvread("dataQ5T4TS2.csv");
D25 = csvread("dataQ5T5TS2.csv");

dt = [2e-2 1e-2, 5e-3, 2.5e-3, 1.25e-3];
dt0 = dt./10;
dt2 = [0.2, 0.1, 0.05, 0.025, 0.0125];

C1Ref = @(x, t) ((4.*pi.*0.01.*t).^(-0.5)).*(exp((-(x - t).^2)/(4.*0.01.*t)));

E = zeros(3,5);
E(1,1) = sum((D01(end,:) - C1Ref(D01(1,:),40)).^2)/1024;
E(1,2) = sum((D02(end,:) - C1Ref(D02(1,:),40)).^2)/1024;
E(1,3) = sum((D03(end,:) - C1Ref(D03(1,:),40)).^2)/1024;
E(1,4) = sum((D04(end,:) - C1Ref(D04(1,:),40)).^2)/1024;
E(1,5) = sum((D05(end,:) - C1Ref(D05(1,:),40)).^2)/1024;

E(2,1) = sum((D11(end,:) - C1Ref(D11(1,:),40)).^2)/1024;
E(2,2) = sum((D12(end,:) - C1Ref(D12(1,:),40)).^2)/1024;
E(2,3) = sum((D13(end,:) - C1Ref(D13(1,:),40)).^2)/1024;
E(2,4) = sum((D14(end,:) - C1Ref(D14(1,:),40)).^2)/1024;
E(2,5) = sum((D15(end,:) - C1Ref(D15(1,:),40)).^2)/1024;

E(3,1) = sum((D21(end,:) - C1Ref(D21(1,:),40)).^2)/1024;
E(3,2) = sum((D22(end,:) - C1Ref(D22(1,:),40)).^2)/1024;
E(3,3) = sum((D23(end,:) - C1Ref(D23(1,:),40)).^2)/1024;
E(3,4) = sum((D24(end,:) - C1Ref(D24(1,:),40)).^2)/1024;
E(3,5) = sum((D25(end,:) - C1Ref(D25(1,:),40)).^2)/1024;

loglog(1./dt0, sqrt(E(1,:)),'r-','linewidth',2);hold on;
loglog(1./dt, sqrt(E(2,:)),'g-','linewidth',2);
loglog(1./dt2, sqrt(E(3,:)),'b-','linewidth',2);
loglog(1./dt0, sqrt(E(1,:)),'r.','markersize',20);
loglog(1./dt, sqrt(E(2,:)),'g.','markersize',20);
loglog(1./dt2, sqrt(E(3,:)),'b.','markersize',20);
loglog([200 2000],5./[200 2000],'k-');
loglog(1./dt2(2:4),2*(dt2(2:4).^2),'k-');
text(500, 1.5e-2, "Slope = -1",'rotation',-46);
text(18, 1e-2, "Slope = -2",'rotation',-64);
xlabel("1/\Delta t");
ylabel("L_2 Error");
title("Temporal Convergence");
legend(["Explicit Euler","Implicit Euler","Crank-Nicolson"],'location','northeast')
fontsize(gca,16,"points");
xlim([3e0 3e4])

%% Problem 6
clc;clear all;close all;
R = zeros(3, 5);

N = [256, 512, 1024, 2048, 4096];
R(1,:) = [1.26e-2, 2.44e-2, 4.76e-2, 9.41e-2, 0.186];
R(2,:) = [8.776e-2, 0.175, 0.351, 0.711, 1.424];
R(3,:) = [9.05e-2, 0.181, 0.370, 0.739, 1.486];

figure('position',[50 50 1000 400])
subplot(1,2,1)
loglog(N, R(1,:),'r-','linewidth',2);hold on;
loglog(N, R(2,:),'g-','linewidth',2);
loglog(N, R(3,:),'b-','linewidth',2);
loglog(N, R(1,:),'r.','markersize',20);
loglog(N, R(2,:),'g.','markersize',20);
loglog(N, R(3,:),'b.','markersize',20);
loglog(N(2:4), 1e-4*N(2:4),'k-')
text(0.75*(2^10), 1e-1, "slope = 1", 'rotation',40)
xlabel("N Points");
title("Time to solution for fixed \Delta t")
ylabel("Time to Solution");
xticks(N)
xlim([2^7, 2^13])
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('2^{%d}', log2(x)), xticks, 'UniformOutput', false));
fontsize(gca,16,"points");

R = R./15000;
subplot(1,2,2);
loglog(N, R(1,:),'r-','linewidth',2);hold on;
loglog(N, R(2,:),'g-','linewidth',2);
loglog(N, R(3,:),'b-','linewidth',2);
loglog(N, R(1,:),'r.','markersize',20);
loglog(N, R(2,:),'g.','markersize',20);
loglog(N, R(3,:),'b.','markersize',20);
loglog(N(2:4), 1e-8*N(2:4),'k-')
text(0.75*(2^10), 1e-5, "slope = 1", 'rotation',40)
legend("Explicit Euler","Implicit Euler","Crank-Nicolson","location","southeast")
xlabel("N Points");
ylabel("Time per time-step");
title("Time per time-step")
xticks(N)
xlim([2^7, 2^13])
set(gca, 'XTickLabel', arrayfun(@(x) sprintf('2^{%d}', log2(x)), xticks, 'UniformOutput', false));
fontsize(gca,16,"points");