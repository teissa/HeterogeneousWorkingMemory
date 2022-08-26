clear
close all

% parameterize the lerning function for the potential (directly computing potential not using Pest)
Nsim =3;   % number of sims per amplitude
kap =5;  % amplitude of von mises to be used
sh = 1/4;   % shift
sc = 5;     % scaling
% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
colormap(cmap)
colormap(cmap)
dx = xsamp(2)-xsamp(1);
% use the random function to sample from the above distribution and plot
randind=[148, 52, 250];
randinp = xsamp(randind); % this identifies the angle theta for each simulation
xfs = linspace(-pi,pi,256); %angles discretized to find distortion

%%
close all

D=0.05;
pot = 1/2/pi+0*xsamp;
pt=pot;
prob=0*xsamp;
f=prob;

for k=1:Nsim
        inp = randinp(k);

    figure; hold on; viscircles([0,0],1,'Color','k','LineWidth',5); 
    for j=1:k
        if j<k
            plot(cos(randinp(j)),sin(randinp(j)),'o','MarkerEdgeColor',cmap(randind(j),:),'MarkerSize',40,'LineWidth',10); 
        else
            plot(cos(randinp(j)),sin(randinp(j)),'o','MarkerEdgeColor',cmap(randind(j),:),'MarkerFaceColor',cmap(randind(j),:),'MarkerSize',40,'LineWidth',10); 
        end
    end
        xticks([]); yticks([]);

    f=f+sc*(sh+exp(kap*(cos(xsamp-inp)-1)))/k;
    prob=f/sum(f);
%     prob=prob+p;
%     prob=prob/k;
      figure, hold on,  xlim([-pi pi]); yticks([]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',40,'interpreter','latex');
ylabel('$P(\theta_j|\theta_{1:j-1})$','fontsize',30,'interpreter','latex');
plot(xsamp,prob,'linewidth',5,'Color','k')
%   for j=1:k
%       if j<k
%         plot(randinp(j),0.015,'o','MarkerEdgeColor',cmap(randind(j),:),'MarkerSize',40);
%       else
%          plot(randinp(j),0.015,'o','MarkerEdgeColor',cmap(randind(j),:),'MarkerFaceColor',cmap(randind(j),:),'MarkerSize',40,'LineWidth',20); 
%       end
% 
%   end  
    
figure, hold on,  xlim([-pi pi]); yticks([]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',40,'interpreter','latex');
ylabel('$U(\theta|\theta_{1:j-1})$','fontsize',30,'interpreter','latex');


    pt = pt+sc*(sh-exp(kap*(cos(xsamp-inp)-1)));
    pot = pt/sum(pt);
    plot(xsamp,pot,'linewidth',5,'Color','k');
%     for j=1:k
%         plot(randinp(j),1,'.','Color',cmap(randind(j),:),'MarkerSize',100);
%     end
end
set(gca,'FontSize',40); set(gca,'TickLabelInterpreter','Latex')

%%
n=4;
A=1;
pr=@(x)A/n*cos(n*x);
ut=@(x)-A/n*cos(n*x);

figure; plot(xsamp,pr(xsamp),'k','LineWidth',5); set(gca,'box','off')
 xlim([-pi pi]); yticks([]); ylim([0 0.015]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$P(\theta_j|\theta_{1:j-1})$','fontsize',30,'interpreter','latex');

figure; plot(xsamp,ut(xsamp),'k','LineWidth',5); set(gca,'box','off')
 xlim([-pi pi]); yticks([]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('$U(\theta_j|\theta_{1:j-1})$','fontsize',30,'interpreter','latex');

 figure; hold on; viscircles([0,0],1,'Color','k','LineWidth',5); 
plot(cos(0),sin(0),'o','MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',40,'LineWidth',10); 

        xticks([]); yticks([]);
