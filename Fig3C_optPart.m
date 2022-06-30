clear
close all

load hetP_varDT

for k=1:3
figure; hold on; 
pcolor(As,ns,rdist{2,k}); shading flat; colormap(flipud(hot)); caxis([0 0.6]); colorbar
 [~,mind]=min(rdist{2,k});
 plot(As,mind,'m.','MarkerSize',20)
axis('tight'); xlim([0.1 2])
 title(['$T=$' num2str(Ts(k))],'Interpreter','Latex')
 set(gca,'fontsize',24);
set(gca,'ticklabelinterpreter','Latex')
xlabel('Amplitude','fontsize',30,'interpreter','latex');
ylabel('Wells','fontsize',30,'interpreter','latex');
end