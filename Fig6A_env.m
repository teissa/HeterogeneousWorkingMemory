close all
clear all
% now build a discretized version of the prior and check sampling from it
sampres = 360;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values

psamp=0.1*ones(1,sampres); psamp(1:10)=1; psamp(80:100)=1; psamp(170:190)=1; psamp(260:280)=1; psamp(350:360)=1;
psamp(1)=0; psamp(end)=0;
figure; hold on; 
plot(xsamp,psamp)
cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
%   cmap=colormap('hsv');
%    cmap=cmap([226:256,1:225],:);
colormap(cmap)
c =1:sampres;
fill(xsamp,psamp,c)
ylim([0 1.2]);  yticks=([]);
xlim([-pi pi]); xticks([-pi/2,0,pi/2]); xticklabels({'-90','0','90'})
set(gca,'fontsize',24);set(gca, 'TickLabelInterpreter','Latex');
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$P_{\rm env}(\theta)$','fontsize',30,'interpreter','latex');


