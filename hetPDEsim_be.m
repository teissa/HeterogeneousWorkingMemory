function [bpost,mse,prb] = hetPDEsim_be(s,A,n,dt,Tim,N, tar,Ntrials, resp)
% simulate the evolution of the probability distribution associated with
% the position of a particle evolving in a potential on a periodic domain
% from -pi to pi. here we use a periodic potential
% define model parameters s = standard deviation of noise
% A = amplitude of potential, n = frequency of potential, dt= timestep, T=
% delays Don't make dt too big! N= number of spatial steps (don't make big)
%tar= targets, Ntrials=number of trials
% need to worry about CFL condition : dt/dx must be sufficiently small
dx = 2*pi/N;    % thus create space step
x = linspace(-pi,pi-dx,N);  % and spatial mesh

U = -(A/n)*cos(n*x);    % formulate discretely meshed potential with params
Uf=[U(2:end),U(1)]; Ub=[U(end), U(1:end-1)];
Uprox=(Uf-Ub)/2/dx;   % use centered difference to approx derivative (2nd order)

% now create the diffusion matrix for simulation
omain = ones(N,1); oof = ones(N-1,1);   % first diffusion
% D = -s^2*diag(omain)+s^2*diag(oof,1)/2+s^2*diag(oof,-1)/2;
% D(1,N)=s^2/2;   D(N,1)=s^2/2;   % corner entries for periodicity
D = -s*diag(omain)+s*diag(oof,1)/2+s*diag(oof,-1)/2;
D(1,N)=s/2;   D(N,1)=s/2;   % corner entries for periodicity
D = D/dx/dx;

I = eye(N); % identity matrix

% make an advection matrix
Ad = diag(Uprox(2:end),1)/2/dx-diag(Uprox(1:end-1),-1)/2/dx;
Ad(1,N) = -Uprox(end)/2/dx; Ad(N,1) = Uprox(1)/2/dx;

% now make the updating backward euler matrix to divide by
M = I-dt*D-dt*Ad; iM = inv(M);


for sm=1:Ntrials %compute for each trial

    t=0;  
    p = zeros(N,1); 
    [~,strt] = min(abs(tar(sm)-x));
    p(strt)=1/dx;  % initialize delta mass at target;
    T=Tim(sm);
    while t<T
        % construct new advection matrix first
        %p = iM*p;
        p = max(iM*p,1e-15); p=p/sum(p)/dx;
        t=t+dt;
    end
    
    prb(sm,:)=p*dx;

    % compute probability of a particular response given this model (for use in
    % Bayesian model/parameter inference)
    xresp = resp(sm);   % subject response
    [~,xind] = min(abs(xresp-x));
    bpost(sm) = prb(sm,xind);
    
    % alternatively compute the mean squared error
    mse(sm)=dx*sum(p'.*(x-xresp).^2);
end
end