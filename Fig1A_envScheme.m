clear all
close all

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 4;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 1; % amplitude of the prior
n=4;
A=1;
pr = @(x) exp(Am*cos(m*x)); %prior function
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
Nsim = 1e5;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values

figure(1), %histogram(randinp, 'Normalization', 'pdf'); %plot distribution of samples
hold on, plot(xsamp,pr(xsamp)/nf,'k','LineWidth',2);
prb=pr(xsamp)/nf;prb(1)=0;prb(end)=0;
cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
%   cmap=colormap('hsv');
%    cmap=cmap([226:256,1:225],:);
colormap(cmap)
c =1:sampres;
fill(xsamp,prb,c)
% and compare to true distribution
% plot([0 0],[0 0.5],'c','LineWidth',2); plot([-pi/2 -pi/2],[0 0.5],'c','LineWidth',2);
% plot([pi/2 pi/2],[0 0.5],'c','LineWidth',2); plot([-pi -pi],[0 0.5],'c','LineWidth',2);
% plot([pi pi],[0 0.5],'c','LineWidth',2); 
ylim([0 0.5])
xlim([-pi pi]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$P_{\rm env}(\theta)$','fontsize',30,'interpreter','latex');

%%
figure; hold on;
r=1; x=0; y=0;
th = -pi:pi/50:pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
h = plot(xunit, yunit);
cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
%   cmap=colormap('hsv');
%    cmap=cmap([226:256,1:225],:);
%cmap=linspecer(360,'sequential');
colormap(cmap)
c=1:length(xunit);
fill(xunit,yunit,c)

