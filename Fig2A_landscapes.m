clear all
close all

n=0;
A=1;
ut=@(x)-A*cos(n*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U_0(\theta)$','fontsize',30,'interpreter','latex');

%%
%close all

n=4;
A=1;
ut=@(x)-A/n*cos(n*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U_4(\theta)$','fontsize',30,'interpreter','latex');

%%
%close all

n=8;
A=1;
ut=@(x)-A/n*cos(n*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U_8(\theta)$','fontsize',30,'interpreter','latex');
