%%  neural field weight examples

clear
close all

% initialize model and simulation parameters
Nx = 1000;   % grid points
Nti = 500;   % input time points
Nt = 10000;   % delay time points
L = 18;      % half length of the grid
Lsc=10; % scaling of half length
Ti = 50;       % input time
I0 = 1;       % Gaussian input
xi = 0;        % initial input location
si = 1;        % spread of Gaussian input
A = 0.35;       % inhibitory strength
sw = 3;        % inhibitory width
T = 1000;      % run time
h = 0.1;     % threshold
s = 0.4;        % strength of weight heterogeneity
n =40;  % frequency of heterogeneity
ep = 0.5;
dx = 2*L/Nx; x = linspace(-L,L,Nx+1); xsc= x.*Lsc;
dt = T/Nt; timey = linspace(0,(T+Ti)/100,Nti+Nt+1);

% setup a matrix which is the weight function (with a periodic
% heterogeneity)
[X,Y] = meshgrid(x,x); 
%gaussian
W_het= dx*(1+s*cos(pi*n*Y/180)).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));
W_het_E= dx*(1+s*cos(pi*n*Y/180)).*(exp(-(X-Y).^2));
W_het_I= dx*(1+s*cos(pi*n*Y/180)).*(-A*exp(-(X-Y).^2/sw^2));
W_hom= dx*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));
W_hom_E= dx*(exp(-(X-Y).^2));
W_hom_I= dx*(-A*exp(-(X-Y).^2/sw^2));

xsamp=xsc*(pi/180);
figure; hold on;
for k=1:20:length(xsamp)
    plot(xsamp,W_het_E(k,:),'color',[1,0.75,0.75],'linewidth',1)
    plot(xsamp,W_het_I(k,:),'color',[0.75,0.75,1],'linewidth',1)
end
plot(xsamp,W_het_E(500,:),'color','r','linewidth',5)
plot(xsamp,W_het_I(500,:),'color','b','linewidth',5)
xlim([-pi/2 pi/2]); ylim([-0.03 0.06]); yticks([]);
xticks([-pi/4,0, pi/4]); xticklabels({'-45','0','45'});
xlabel('Preferred $\theta$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');


figure; hold on;
for k=1:20:length(xsamp)
    plot(xsamp,W_hom_E(k,:),'color',[1,0.75,0.75],'linewidth',1)
    plot(xsamp,W_hom_I(k,:),'color',[0.75,0.75,1],'linewidth',1)
end
plot(xsamp,W_hom_E(500,:),'color','r','linewidth',5)
plot(xsamp,W_hom_I(500,:),'color','b','linewidth',5)
xlim([-pi/2 pi/2]); ylim([-0.03 0.06]); yticks([])
xticks([-pi/4,0, pi/4]); xticklabels({'-45','0','45'});

xlabel('Preferred $\theta$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');
