clear
close all


%% Example of actibvation of the initial trial

% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values

cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
colormap(cmap)

figure; hold on;
for k=1:10:size(cmap,1)
    plot(xsamp(k),0,'.','MarkerSize',40,'Color',cmap(k,:))
end

% initialize model and simulation parameters
Nx = 1000;   % grid points
Nti = 500;   % input time points
Nt = 10000;   % delay time points
L = 18;      % half length of the grid
Lsc=10; % scaling of half length
Ti = 50;       % input time
I0 = 0.3;       % Gaussian input
si = 1;        % spread of Gaussian input
xi=0; %input location
A = 0.35;       % inhibitory strength
sw = 3;        % inhibitory width
T = 1000;      % run time
h = 0.1;     % threshold
s = 0.5;        % strength of weight heterogeneity
kap = 0.3;  % scaling of bump profile for learning
n =40;  % frequency of heterogeneity
ep = 0.4;
dx = 2*L/Nx; x = linspace(-L,L,Nx+1); xsc= x.*Lsc;
dt = T/Nt; timey = linspace(0,(T+Ti)/100,Nti+Nt+1);

% setup a matrix which is the weight function (with a periodic
% heterogeneity)
[X,Y] = meshgrid(x,x); 
het=0*Y;

%gaussian
W_learn= dx*(1+s*het).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));

% initialize the neural field activity variable
U_learn = zeros(Nx+1,Nti+Nt+1);

% this is the noise filter
F = sqrt(dx)*sqrt(ep)*exp(-abs(X-Y));
cm_learn(1)=xi;

% start by applying an input for 500 ms
I = I0*exp(-(xsc-xi).^2/2/si^2); I=I';
for j=1:Nti
    % update neural field
    Hu_learn = heaviside(U_learn(:,j)-h);   % sample the current active region

    nos = F*randn(Nx+1,1);      % make noise for this timestep

    % complete the current Euler-Maruyama step
    U_learn(:,j+1) = U_learn(:,j)+dt*(-U_learn(:,j)+W_learn*Hu_learn+I)+sqrt(dt)*U_learn(:,j).*nos;
     [junk,mi] = max(U_learn(:,j+1));
    cm_learn(j+1) = xsc(mi);
    
end

%identify bump average
bp=mean(U_learn(:,1:Nti)');

%update learning
 het = 0.99*het + kap*bp; 
W_learn= dx*(1+s*het).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));

plot(xsc*(pi/180),U_learn(:,Nti),'color','k','linewidth',2)
xlim([-pi pi]); %ylim([-0.01 0.04]); 
yticks([]); xticks([-pi/2,0,pi/2]); xticklabels({'-90','0','90'});
xlabel('Preferred $\theta$','fontsize',30,'interpreter','latex');
ylabel('Neural Activity','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');


figure; pcolor(xsc,xsc,W_learn); shading flat
xlabel(' Preferred $\theta$','fontsize',30,'interpreter','latex'); ylabel(' Preferred $\theta$','fontsize',30,'interpreter','latex');
caxis([-0.01 0.03]); colorbar; colormap('hot')
xticks([-90,0,90]); yticks([-90,0,90]); 
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');
%% Example of the secondary trial

xi = (pi/4)*180/pi;% initialize the neural field activity variable
U_learn = zeros(Nx+1,Nti+Nt+1);

% this is the noise filter
F = sqrt(dx)*sqrt(ep)*exp(-abs(X-Y));
cm_learn(1)=xi;

% start by applying an input for 500 ms
I = I0*exp(-(xsc-xi).^2/2/si^2); I=I';
for j=1:Nti
    % update neural field
    Hu_learn = heaviside(U_learn(:,j)-h);   % sample the current active region

    nos = F*randn(Nx+1,1);      % make noise for this timestep

    % complete the current Euler-Maruyama step
    U_learn(:,j+1) = U_learn(:,j)+dt*(-U_learn(:,j)+W_learn*Hu_learn+I)+sqrt(dt)*U_learn(:,j).*nos;
     [junk,mi] = max(U_learn(:,j+1));
    cm_learn(j+1) = xsc(mi);
    
end

%identify bump average
bp=mean(U_learn(:,1:Nti)');

%update learning
 het = 0.99*het + kap*bp; 
W_learn= dx*(1+s*het).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));

figure; hold on;
for k=1:10:size(cmap,1)
    plot(xsamp(k),0,'.','MarkerSize',40,'Color',cmap(k,:))
end
plot(xsc*(pi/180),U_learn(:,Nti),'color','k','linewidth',2)
xlim([-pi pi]); %ylim([-0.01 0.04]); 
yticks([]); xticks([-pi/2,0,pi/2]); xticklabels({'-90','0','90'});
xlabel('Preferred $\theta$','fontsize',30,'interpreter','latex');
ylabel('Neural Activity','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');


figure; pcolor(xsc,xsc,W_learn); shading flat
xlabel(' Preferred $\theta$','fontsize',30,'interpreter','latex'); ylabel(' Preferred $\theta$','fontsize',30,'interpreter','latex');
caxis([-0.01 0.03]); colorbar
xticks([-90,0,90]); yticks([-90,0,90]); 
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');

%% example weights for the static heterogeneous case

% initialize model and simulation parameters
Nx = 1000;   % grid points
Nti = 500;   % input time points
Nt = 10000;   % delay time points
L = 18;      % half length of the grid
Lsc=10; % scaling of half length
Ti = 50;       % input time
I0 = 1;       % Gaussian input
xi = (-pi/2)*180/pi;        % initial input location
si = 1;        % spread of Gaussian input
A = 0.4;       % inhibitory strength
sw = 3;        % inhibitory width
T = 1000;      % run time
h = 0.1;     % threshold
s = 0.4;        % strength of weight heterogeneity
n =40;  % frequency of heterogeneity
ep = 0.5;
dx = 2*L/Nx; x = linspace(-L,L,Nx+1); xsc= x.*Lsc;
dt = T/Nt; timey = linspace(0,(T+Ti)/100,Nti+Nt+1);

% initialize the neural field activity variable
U_het = zeros(Nx+1,Nti+Nt+1);

% setup a matrix which is the weight function (with a periodic
% heterogeneity)
[X,Y] = meshgrid(x,x); 
%gaussian
W_het= dx*(1+s*cos(pi*n*Y/180)).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));

% this is the noise filter
F = sqrt(dx)*sqrt(ep)*exp(-abs(X-Y));
cm_het(1)=xi;
% start by applying an input for 500 ms
I = I0*exp(-(xsc-xi).^2/2/si^2); I=I';
for j=1:Nti
    % update neural field
    Hu_het = heaviside(U_het(:,j)-h);   % sample the current active region

    nos = F*randn(Nx+1,1);      % make noise for this timestep

    % complete the current Euler-Maruyama step
    U_het(:,j+1) = U_het(:,j)+dt*(-U_het(:,j)+W_het*Hu_het+I)+sqrt(dt)*U_het(:,j).*nos;
     [junk,mi] = max(U_het(:,j+1));
    cm_het(j+1) = xsc(mi);

end



figure; pcolor(xsc,xsc,W_het); shading flat
xlabel(' Preferred $\theta$','fontsize',30,'interpreter','latex'); ylabel(' Preferred $\theta$','fontsize',30,'interpreter','latex');
caxis([-0.01 0.03]); colorbar
xticks([-90,0,90]); yticks([-90,0,90]); 
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');

%% connectivity in the static heterogeneous case

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


xsamp=xsc*(pi/180);
figure; hold on;
for k=1:10:length(xsamp)
    plot(xsamp,W_het(k,:),'color','k','linewidth',1)
%     plot(xsamp,W_het_E(k,:),'color',[1,0.75,0.75],'linewidth',1)
%     plot(xsamp,W_het_I(k,:),'color',[0.75,0.75,1],'linewidth',1)
end
xlim([-pi pi]); ylim([-0.03 0.06]); yticks([]);
xticks([-pi/2,0, pi/2]); xticklabels({'-90','0','90'});
xlabel('Preferred $\theta$','fontsize',30,'interpreter','latex');
ylabel('$W$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');



% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values

cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
colormap(cmap)

for k=1:10:size(cmap,1)
    plot(xsamp(k),0,'.','MarkerSize',40,'Color',cmap(k,:))
end


