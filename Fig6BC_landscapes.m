%%
clear all
close all

%%

n1=0;
A1=1;
ut=@(x)-A1*cos(n1*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');


%%

n1=4;
A1=1;
ut=@(x)-A1*cos(n1*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');


%%
n1=4;
A1=1;
A2=0.5;
n2=8;
ut=@(x)-A1*cos(n1*x)-A2*cos(n2*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');

%%
n1=4;
A1=1;
offs=pi/2;
ut=@(x)-A1*cos(n1*x-offs);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
figure;
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');

%%



n1=4;
A1=1;
ut=@(x)-A1*cos(n1*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
dx=1/sampres;
sh=0.25;
sc=1;
kap=8;
inp=-4*pi/6;

learnPot=zeros(1,sampres);
learnPot = learnPot+(sc*(sh-exp(kap*(cos(xsamp-inp)-1))));
learnPot = learnPot/sum(learnPot);

figure;
hold on, plot(xsamp,learnPot,'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');



learnPot=ut(xsamp);
learnPot = learnPot+(sc*(sh-exp(kap*(cos(xsamp-inp)-1))));
learnPot = learnPot/sum(learnPot);

figure;
hold on, plot(xsamp,learnPot,'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');


%%



n1=4;
A1=1;
ut=@(x)-A1*cos(n1*x);
% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
dx=1/sampres;
sh=0.25;
sc=1;
kap=8;
inp=-4*pi/6;
dis=pi/3;

learnPot=zeros(1,sampres);
learnPot = learnPot+(sc*(sh-exp(kap*(cos(xsamp-inp)-1))));
learnPot = learnPot+(sc*(sh-exp(kap*(cos(xsamp-dis)-1))));
learnPot = learnPot/sum(learnPot);


figure;
hold on, plot(xsamp,learnPot,'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');



learnPot=ut(xsamp);
learnPot = learnPot+(sc*(sh-exp(kap*(cos(xsamp-inp)-1))));
learnPot = learnPot+(sc*(sh-exp(kap*(cos(xsamp-dis)-1))));
learnPot = learnPot/sum(learnPot);

figure;
hold on, plot(xsamp,learnPot,'k','LineWidth',2);  yticks([])
xlim([-pi pi])% and compare to true distribution
xticks([-pi/2, 0, pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target($\theta$)','fontsize',30,'interpreter','latex');
ylabel('$U(\theta)$','fontsize',30,'interpreter','latex');
