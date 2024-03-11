clear all; clc;
f = figure;
f.Position = [488,242,10000,350];

subplot(1,4,1);
n = 0.119:0.019:2;
data = load('Data\noise_gaussian_diffdelta_N2.mat');
plot(n,data.obj2 - data.obj1,'-.',"Color",'r','LineWidth',1.1); hold on;
plot(n,data.obj3 - data.obj2,"Color",'b','LineWidth',1.1); hold on;
plot(n,data.obj4 - data.obj3,'--',"Color",'g','LineWidth',1.1); hold on;
lgd=legend({'$\mathcal{J}_{\mathrm{max}}^{(ii)}-\mathcal{J}_{\mathrm{max}}^{(i)}$', ...
    '$\mathcal{J}_{\mathrm{max}}^{(iii)}-\mathcal{J}_{\mathrm{max}}^{(ii)}$', ...
    '$\mathcal{J}_{\mathrm{max}}^{(iv)}-\mathcal{J}_{\mathrm{max}}^{(iii)}$'},'Interpreter','latex','Location','northwest');
set(gca,'XLim',[0.1 2]);
fontsize(lgd,12,'points');
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('$\Delta$','Interpreter','latex','FontSize',20);
grid minor;

subplot(1,4,2);
n = -2.94:0.06:3;
data = load('Data\noise_gaussian_diffmean_N2.mat');
plot(n,data.obj2 - data.obj1,'-.',"Color",'r','LineWidth',1.1); hold on;
plot(n,data.obj3 - data.obj2,"Color",'b','LineWidth',1.1); hold on;
plot(n,data.obj4 - data.obj3,'--',"Color",'g','LineWidth',1.1); hold on;
set(gca,'XLim',[-3 3]);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('$\mu$','Interpreter','latex','FontSize',20);
grid minor;

subplot(1,4,3);
n = 0.01:0.01:1;
data = load('Data\noise_GM_diffw_N2.mat');
plot(n,data.obj2 - data.obj1,'-.',"Color",'r','LineWidth',1.1); hold on;
plot(n,data.obj3 - data.obj2,"Color",'b','LineWidth',1.1); hold on;
plot(n,data.obj4 - data.obj3,'--',"Color",'g','LineWidth',1.1); hold on;
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('$w$','Interpreter','latex','FontSize',20);
grid minor;

subplot(1,4,4);
n = 0.218:0.018:2;
data = load('Data\noise_beta_diffa_N2.mat');
plot(n,data.obj2 - data.obj1,'-.',"Color",'r','LineWidth',1.1); hold on;
plot(n,data.obj3 - data.obj2,"Color",'b','LineWidth',1.1); hold on;
plot(n,data.obj4 - data.obj3,'--',"Color",'g','LineWidth',1.1); hold on;
set(gca,'XLim',[0.2 2]);
ax = gca;
ax.XAxis.FontSize = 12;
ax.YAxis.FontSize = 12;
xlabel('$a$','Interpreter','latex','FontSize',20);
grid minor;

saveas(gcf,'Plot.svg');