clear
close all

%create the prior and parameterize it to build the sampling 

D = 0;   % noise corresponding to corrosion of parametric WM rep, diffusion
T = 5;      % delay time (seconds)
dt = 0.01;  % timestep of stochastic sims
nt = round(T/dt)+1; % number of timesteps per sim

% now parameterize the prior
m = 4;  % number of peaks on the periodic prior
Am = 2; % amplitude of the prior
pr = @(x) exp(Am*cos(m*x)); %prior function
nf = integral(pr,-pi,pi); %integral of the prior (for normalizing)
cpr = @(x)integral(pr,-pi,x)/nf; %cumulative probability distribution function

% build sampling
Nsim = 500;   % number of sims per amplitude

% now build a discretized version of the prior and check sampling from it
sampres = 5e3;  % discretization of the prior function
xsamp = linspace(-pi,pi,sampres); %angle values
Fsamp = zeros(1,sampres);
for j=1:sampres, Fsamp(j) = cpr(xsamp(j)); end %builds the cumulative probability 

% now we make an object that we can sample from 
pdist = makedist('PiecewiseLinear','x',xsamp,'Fx',Fsamp);

% use the random function to sample from the above distribution and plot
randinp = random(pdist,1,Nsim); % this identifies the angle theta for each simulation

kap =8;  % amplitude of von mises to be used
sh = 1/4;   % shift
sc = 5;     % scaling
xs = linspace(-pi,pi,1000); dx = xs(2)-xs(1);

for rr=1:10
    r=randperm(Nsim);


    %initialize parameters
    pot = 1/2/pi+0*xs; %initial potential landscape (flat)
    rout = zeros(Nsim,1); %responses at end of trial
    rdist = 0; %distortion
    for l=1:Nsim
        inp = randinp(r(l));
        x=inp;
        pvar(l,:)=pot;

        pgrad = gradient(pot,dx); %find the derivative of the current potential landscape

        % simulate the delayed response trial
        for j=1:nt-1, x=mod(x-dt*interp1(xs,pgrad,x)+sqrt(dt*2*D)*randn+pi,2*pi)-pi; end
        rout(l) = x;%response at end with drift and diffusion

        % add the distortion value for this trial
        rdist = (rdist+(1-cos(x-inp)));  % running average of distortion
        trial_dist(l)=(1-cos(x-inp));
        mdist(l)=rdist/l;

        %update the potential landscape using the previous trial input
        pot = pot+sc*(sh-exp(kap*(cos(xs-inp)-1)))/l;
        pot = pot/dx/sum(pot);

    end

% plot the running distortion average
figure(1), hold on ; xlim([1 Nsim]); ylim([0 0.5])
plot(mdist,'LineWidth',5);
xlabel('Trial','interpreter','latex','fontsize',30);
ylabel('Distortion','interpreter','latex','fontsize',30);
set(gca,'fontsize',30); set(gca,'ticklabelinterpreter','Latex')


figure(2), hold on,  xlim([-pi pi]); yticks([])
set(gca,'fontsize',24); set(gca,'ticklabelinterpreter','Latex')
xlabel('$\theta$','fontsize',30,'interpreter','latex');
ylabel('Potential','fontsize',30,'interpreter','latex');
plot(xs,pot,'linewidth',4);
set(gca,'FontSize',30); set(gca,'TickLabelInterpreter','Latex')
end
