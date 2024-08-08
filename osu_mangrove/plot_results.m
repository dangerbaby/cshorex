

close all
fs = 12;
iprint = 0;
figure; 
plot(out.x,out.Hrms,'linewidth',2)
ylabel('$H_{rms}[m]$','fontsize',fs,'interpreter','latex')
xlabel('$x[m]$','fontsize',fs,'interpreter','latex')
title(['Wave Height: ',strrep(in.name,'_','-')],'fontsize',fs,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

figure;
plot(out.x,out.eta,'linewidth',2);hold all

fill([out.x out.x(end) out.x(1)],[in.zb in.zb(1) in.zb(1)],[.8 .8 .5])
fill([out.x fliplr(out.x)],[out.eta fliplr(in.zb)],[.7 .7 1])
% fill([out.x(find(out.iswash==0)) fliplr(out.x(find(out.iswash==0)))],...
%      [in.zb(find(out.iswash==0)) fliplr(out.eta(find(out.iswash==0)))],[.7 .7 1],'facealpha',.6);
% fill([out.x(find(out.iswash==1)) fliplr(out.x(find(out.iswash==1)))],...
%      [in.zb(find(out.iswash==1)) fliplr(out.eta(find(out.iswash==1)))],[.7 .7 1],'facealpha',.2);

plot(out.x,0*out.x+in.swlbc,'linewidth',1);hold all
plot(out.runup_2p_x,out.runup_2p,'rs','markersize',10,'markerfacecolor','r')  
title(['Mean Free Surface: ',strrep(in.name,'_','-')],'fontsize',fs,'interpreter','latex')
ylabel('$\eta [m]$','fontname','times','fontsize',fs,'interpreter','latex')
xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
plot(in.x,in.zb,'k')
ylim([-1 .5])
%axis([0 150 -1 2])
set(gca,'TickLabelInterpreter','latex')
set(gcf,'Position',[672 545 800 400])
if iprint;print -dpng eta.png;end  

figure;
plot(out.x,out.Ef);hold all
ylabel('$E_f$','fontsize',fs,'interpreter','latex')
xlabel('$x[m]$','fontsize',fs,'interpreter','latex')
set(gca,'TickLabelInterpreter','latex')

if max(max(abs(out.V)))>.001
    figure;
  plot(out.x,out.V);hold all
  ylabel('$V[m/s]$','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontsize',fs,'interpreter','latex')
  set(gca,'TickLabelInterpreter','latex')
end

if in.iprofl 
  figure;
  plot(out.x,out.zb);hold all
  ylabel('$zb[m]$','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontsize',fs,'interpreter','latex')
  set(gca,'TickLabelInterpreter','latex')
end

if isfield(out.sed,'Vb')
  figure
  plot(out.x,out.sed(end).aveqby,'b','linewidth',2);hold all

  plot(out.x,out.sed(end).aveqsy,'r','linewidth',2);hold all
  plot(out.x,out.sed(end).aveqy,'k','linewidth',2);hold all
  title('Longshore Sediment Transport ',...
        'fontname','times','fontsize',fs,'interpreter','latex')
  ylabel('$Q_y[m^2/s]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  if iprint;print -dpng longshore_tranport.png;end

  figure
  plot(out.x,0*out.sed(end).aveqbx,'k','linewidth',2);hold all
    plot(out.x,out.sed(end).aveqbx,'b','linewidth',2);hold all
  plot(out.x,out.sed(end).aveqbx+out.sed(end).stdqbx,'b--','linewidth',2);hold all
  plot(out.x,out.sed(end).aveqbx-out.sed(end).stdqbx,'b--','linewidth',2);hold all
  plot(out.x,out.sed(end).aveqsx,'r','linewidth',2);hold all
  plot(out.x,out.sed(end).aveqsx+out.sed(end).stdqsx,'r--','linewidth',2);hold all
  plot(out.x,out.sed(end).aveqsx-out.sed(end).stdqsx,'r--','linewidth',2);hold all
  title('Cross Shore Sediment Transport ',...
        'fontname','times','fontsize',fs,'interpreter','latex')
  ylabel('$Q_x[m^2/s]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  if iprint;print -dpng cross_shore_tranport.png;end

  figure
  ind = [200];
  %plot(out.x,out.sed(end).Cb(ind,:),'b','linewidth',2);hold all
  plot(out.sed(end).time(:,ind),out.sed(end).Cb(:,ind),'b','linewidth',2);hold all
    plot(out.sed(end).time(:,ind),out.sed(end).Vbe(:,ind)/.01,'b--','linewidth',2);hold all
  plot(out.sed(end).time(:,ind),out.sed(end).Cs0(:,ind),'r','linewidth',2);hold all
  %plot(out.x,out.sed(end).Cs0(ind,:),'r','linewidth',2);hold all
  title('Sediment Conc ',...
        'fontname','times','fontsize',fs,'interpreter','latex')
  ylabel('$C[m^3/m^3]$','fontname','times','fontsize',fs,'interpreter','latex')
  xlabel('$x[m]$','fontname','times','fontsize',fs,'interpreter','latex')
  if iprint;print -dpng sed_conc.png;end

end


% if in.iprofl 
%   figure;
%   hh(1) = plot(out.x,out.sed.aveqbx);hold all
%   hh(2) = plot(out.x,out.sed.aveqsx);hold all
%   hh(3) = plot(out.x,out.sed.aveqx);hold all
%   legend(hh,'Bedload','Suspended Load','Total','location','northwest')
%   ylabel('$q[m^2/s]$','interpreter','latex')
%   xlabel('$x[m]$','interpreter','latex')
%   set(gca,'TickLabelInterpreter','latex')
% end

