clear
close all

%% Neural Field Bump Example

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
T = 500;      % run time
h = 0.1;     % threshold
s = 0.4;        % strength of weight heterogeneity
n =40;  % frequency of heterogeneity
ep = 0.6;
dx = 2*L/Nx; x = linspace(-L,L,Nx+1); xsc= x.*Lsc;
dt = T/Nt; timey = linspace(0,(T+Ti)/100,Nti+Nt+1);

% initialize the neural field activity variable
U_hom = zeros(Nx+1,Nti+Nt+1);

% setup a matrix which is the weight function (with a periodic
% heterogeneity)
[X,Y] = meshgrid(x,x); 
%gaussian
W_homG= dx*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));


%exponential
W_hom = dx*(1-abs(X-Y)).*exp(-abs(X-Y));

% this is the noise filter
F = sqrt(dx)*sqrt(ep)*exp(-abs(X-Y));
cm_hom(1)=xi;
% start by applying an input for 500 ms
I = I0*exp(-(x-xi).^2/2/si^2); I=I';
for j=1:Nti
    % update neural field
    Hu_hom = heaviside(U_hom(:,j)-h);   % sample the current active region

    nos = F*randn(Nx+1,1);      % make noise for this timestep
    
     % complete the current Euler-Maruyama step
    U_hom(:,j+1) = U_hom(:,j)+dt*(-U_hom(:,j)+W_homG*Hu_hom+I)+sqrt(dt)*U_hom(:,j).*nos;
     [junk,mi] = max(U_hom(:,j+1));
    cm_hom(j+1) = xsc(mi);
end

% follow this with a 5s (5000ms) delay period
for j=1:Nt
    % update neural field
    Hu_hom = heaviside(U_hom(:,Nti+j)-h);   % sample the current active region
    nos = F*randn(Nx+1,1);      % make noise for this timestep

    % complete the current Euler-Maruyama step
    
    U_hom(:,Nti+j+1) = U_hom(:,Nti+j)+dt*(-U_hom(:,Nti+j)+W_homG*Hu_hom)+sqrt(dt)*U_hom(:,Nti+j).*nos;
    [junk,mi] = max(U_hom(:,Nti+j+1));
    cm_hom(Nti+j+1) = xsc(mi);
end

figure; hold on; 
pcolor(xsc,timey,U_hom'); shading flat, colormap(hot); caxis([h 0.6])
plot(cm_hom,timey,'m','linewidth',3)
xlabel('Target','Interpreter','Latex')
ylabel('Time','Interpreter','Latex')
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
axis('tight')
%% Builds Particle Energy Landscape

clear all
close all

D = 0.05;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 1; % amplitude of the prior
n=4;
A=1;
pr = @(x) exp(Am*cos(m*x)); %prior function
ut=@(x)-A/n*cos(n*x);
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function

% build sampling
Nsim = 1e5;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 1e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Nsim); % this identifies the angle theta for each simulation

figure, %histogram(randinp, 'Normalization', 'pdf'); %plot distribution of samples
hold on, plot(xsamp,ut(xsamp),'k','LineWidth',2);  yticks([])
% plot([0 0],[-1 1],'c','LineWidth',2); plot([-pi/2 -pi/2],[-1 1],'c','LineWidth',2);
% plot([pi/2 pi/2],[-1 1],'c','LineWidth',2); plot([-pi -pi],[-1 1],'c','LineWidth',2);
% plot([pi pi],[-1 1],'c','LineWidth',2);
xlim([-pi pi])% and compare to true distribution
xticks([-1,0,1]); xticklabels({'-50','0','50'})
%ylim([0 0.5])
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','latex')
xlabel('Target','fontsize',30,'interpreter','latex');
ylabel('Potential','fontsize',30,'interpreter','latex');
xlim([-1.5 1.5])

%% Schematic of Excitatory- Inhibitory Connectivity

clear
close all


% now build a discretized version of the prior and check sampling from it
sampres = 256;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values

cmap=customcolormap([0 0.2 0.4 0.6 0.8 1],{'#fcbdbd','#fcecbd','#bdfcce','#bdd0fc','#c7bdfc','#fcbdc3'});
colormap(cmap)
xs=[1:16:sampres,sampres];
figure; hold on;
for k=xs
    plot(xsamp(k),0,'.','MarkerSize',70,'Color',cmap(k,:))
end

% %%%%% excitatory connectivity
% th = linspace( 0, pi, 100);
% xsm=xsamp(xs);
% for k=1:4
%     xn=xsm(8-k);
%     r=k/2;
%     x = r*cos(th) + xn;
%     y = r*sin(th);
%     plot(x,y,'k','Linewidth',5-k)
% end
% 
for k=1:30:size(cmap,1)
    plot(xsamp(k)+0.14,-0.001,'.','MarkerSize',40,'Color','k')
end

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
W_hom= dx*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));

%plot(xsc*(pi/180),W_hom(500,:),'color','k','linewidth',2)
xlim([-pi pi]); ylim([-0.005 0.005]); 
yticks([]); xticks([-2,0,2]); xticklabels({'-100','0','100'})
xlabel('Preferred $\theta$','fontsize',30,'interpreter','latex');
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex');
