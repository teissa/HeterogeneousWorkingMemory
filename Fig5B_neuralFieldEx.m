%% examples of the bump in the neural field

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
xi = (0)*180/pi;        % initial input location
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
U_hom=U_het;

% setup a matrix which is the weight function (with a periodic
% heterogeneity)
[X,Y] = meshgrid(x,x); 
%gaussian
W_het= dx*(1+s*cos(pi*n*Y/180)).*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));
W_hom= dx*(exp(-(X-Y).^2)-A*exp(-(X-Y).^2/sw^2));

% this is the noise filter
F = sqrt(dx)*sqrt(ep)*exp(-abs(X-Y));
cm_het(1)=xi;
cm_hom(1)=xi;
% start by applying an input for 500 ms
I = I0*exp(-(xsc-xi).^2/2/si^2); I=I';
for j=1:Nti
    % update neural field
    Hu_het = heaviside(U_het(:,j)-h);   % sample the current active region
    Hu_hom = heaviside(U_hom(:,j)-h);   % sample the current active region

    nos = F*randn(Nx+1,1);      % make noise for this timestep

    % complete the current Euler-Maruyama step
    U_het(:,j+1) = U_het(:,j)+dt*(-U_het(:,j)+W_het*Hu_het+I)+sqrt(dt)*U_het(:,j).*nos;
     [junk,mi] = max(U_het(:,j+1));
    cm_het(j+1) = xsc(mi);
    
        nos = F*randn(Nx+1,1);      % make noise for this timestep

     % complete the current Euler-Maruyama step
    U_hom(:,j+1) = U_hom(:,j)+dt*(-U_hom(:,j)+W_hom*Hu_hom+I)+sqrt(dt)*U_hom(:,j).*nos;
     [junk,mi] = max(U_hom(:,j+1));
    cm_hom(j+1) = xsc(mi);
end

% follow this with a 5s (5000ms) delay period
for j=1:Nt
    % update neural field
    Hu_het = heaviside(U_het(:,Nti+j)-h);   % sample the current active region
    Hu_hom = heaviside(U_hom(:,Nti+j)-h);   % sample the current active region
    nos = F*randn(Nx+1,1);      % make noise for this timestep

    % complete the current Euler-Maruyama step
    U_het(:,Nti+j+1) = U_het(:,Nti+j)+dt*(-U_het(:,Nti+j)+W_het*Hu_het)+sqrt(dt)*U_het(:,Nti+j).*nos;
    [junk,mi] = max(U_het(:,Nti+j+1));
    cm_het(Nti+j+1) = xsc(mi);
    
        nos = F*randn(Nx+1,1);      % make noise for this timestep

    U_hom(:,Nti+j+1) = U_hom(:,Nti+j)+dt*(-U_hom(:,Nti+j)+W_hom*Hu_hom)+sqrt(dt)*U_hom(:,Nti+j).*nos;
    [junk,mi] = max(U_hom(:,Nti+j+1));
    cm_hom(Nti+j+1) = xsc(mi);
end

figure; hold on; 
pcolor(xsc*(pi/180),timey,U_het'); shading flat, colormap(hot); caxis([h 0.6]); axis('tight')
plot(cm_het*(pi/180),timey,'m','linewidth',3); plot([xi*pi/180 xi*pi/180],[0 0.5],'g','linewidth',10);
plot([0 0],[0 timey(end)],'c','linewidth',3); plot([-pi/2 -pi/2],[0 timey(end)],'c','linewidth',3);
plot([pi/2 pi/2],[0 timey(end)],'c','linewidth',3);plot([-pi -pi],[0 timey(end)],'c','linewidth',3); plot([pi pi],[0 timey(end)],'c','linewidth',3)
xticks([-pi/2,0,pi/2]); xticklabels({'-90','0','90'});
xlabel('$\theta$','Interpreter','Latex')
ylabel('Time','Interpreter','Latex')
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','Latex');

figure; hold on; 
pcolor(xsc*(pi/180),timey,U_hom'); shading flat, colormap(hot); caxis([h 0.6]); axis('tight')
plot(cm_hom*(pi/180),timey,'m','linewidth',3); plot([xi*pi/180 xi*pi/180],[0 0.5],'g','linewidth',10);
xticks([-pi/2,0,pi/2]); xticklabels({'-90','0','90'});
xlabel('$\theta$','Interpreter','Latex')
ylabel('Time','Interpreter','Latex')
set(gca,'fontsize',30);set(gca,'ticklabelinterpreter','Latex');
