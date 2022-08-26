clear
close all

load hetP_varDT

% sigma=0.01
for k=1:4
figure; hold on; 
pcolor(As,ns,rdist{1,k}); shading flat; colormap(flipud(hot)); caxis([0 0.6]); colorbar
 [~,mind]=min(rdist{1,k});
 plot(As,mind,'m.','MarkerSize',20)
axis('tight'); xlim([0.1 2])
 title(['$T=$' num2str(Ts(k))],'Interpreter','Latex')
 set(gca,'fontsize',24);
set(gca,'ticklabelinterpreter','Latex')
xlabel('Amplitude','fontsize',30,'interpreter','latex');
ylabel('Troughs','fontsize',30,'interpreter','latex');
end

% sigma=0.05
for k=1:4
figure; hold on; 
pcolor(As,ns,rdist{2,k}); shading flat; colormap(flipud(hot)); caxis([0 0.6]); colorbar
 [~,mind]=min(rdist{2,k});
 plot(As,mind,'m.','MarkerSize',20)
axis('tight'); xlim([0.1 2])
 title(['$T=$' num2str(Ts(k))],'Interpreter','Latex')
 set(gca,'fontsize',24);
set(gca,'ticklabelinterpreter','Latex')
xlabel('Amplitude','fontsize',30,'interpreter','latex');
ylabel('Troughs','fontsize',30,'interpreter','latex');
end

%sigma=0.1
for k=1:4
figure; hold on; 
pcolor(As,ns,rdist{1,k}); shading flat; colormap(flipud(hot)); caxis([0 0.6]); colorbar
 [~,mind]=min(rdist{1,k});
 plot(As,mind,'m.','MarkerSize',20)
axis('tight'); xlim([0.1 2])
 title(['$T=$' num2str(Ts(k))],'Interpreter','Latex')
 set(gca,'fontsize',24);
set(gca,'ticklabelinterpreter','Latex')
xlabel('Amplitude','fontsize',30,'interpreter','latex');
ylabel('Troughs','fontsize',30,'interpreter','latex');
end

for k=1:4
figure; hold on; 
pcolor(As,ns,rdist{3,k}); shading flat; colormap(flipud(hot)); caxis([0 0.6]); colorbar
 [~,mind]=min(rdist{3,k});
 plot(As,mind,'m.','MarkerSize',20)
axis('tight'); xlim([0.1 2])
 title(['$T=$' num2str(Ts(k))],'Interpreter','Latex')
 set(gca,'fontsize',24);
set(gca,'ticklabelinterpreter','Latex')
xlabel('Amplitude','fontsize',30,'interpreter','latex');
ylabel('Troughs','fontsize',30,'interpreter','latex');
end