clear all; clc;

rep_len = 100;
step = (2-0.1)/rep_len;
n = (0.1+step):step:2;

% subplot(1,2,1);
% GHZ
% N = 2
data = load('Data\GHZ_diffdelta_N2.mat');
plot(n,data.objlocal./data.objglobal,'-.',"Color",'r','LineWidth',1.1); hold on;


% N = 3
data = load('Data\GHZ_diffdelta_N3.mat');
plot(n,data.objlocal./data.objglobal,"Color",'b','LineWidth',1.1); hold on;
xlabel('$\Delta$','Interpreter','latex');
ylabel('$\frac{\mathcal{J}(\rho)}{\mathcal{J}_{\mathrm{max}}^{(i)}}$','Interpreter','latex');
legend({'N=2','N=3'});
set(gca,'yminortick','on');
set(gca,'xminortick','on');
set(gca,'XLim',[0.1+step,2]);
saveas(gcf,'Plot.svg');

% This figure is similar to the one above
% subplot(1,2,2);
% % Seqential with or without control
% % N = 2
% data = load('Data\SeqCtrl_diffprior_N2.mat');
% plot(n,data.objlocal./data.objglobal,'-.',"Color",'r','LineWidth',1.1); hold on;
% 
% 
% % N = 3
% data = load('Data\SeqCtrl_diffprior_N3.mat');
% plot(n,data.objlocal./data.objglobal,"Color",'b','LineWidth',1.1); hold on;
% xlabel('$\Delta$','Interpreter','latex');
% ylabel('$\frac{\mathcal{J}_{\mathrm{max}}^{(ii)}}{\mathcal{J}(\rho)}$','Interpreter','latex');
% legend({'N=2','N=3'});