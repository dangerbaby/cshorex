close all

figure;clf;clear hh
dum = 'osv><^dhp';
for i = 1:length(sve)
  hhdum = plot(sve(i).meanu,sve(i).undertow,['r',dum(i)],'markerfacecolor','k','markersize',10);hold on
  hh(i) = hhdum(1);
  hlabs{i} = dnames(i).name;
end
hl = legend(hh,hlabs);
set(hl,'interpreter','latex','fontsize',14,'location','northwest','autoupdate','off')
plot([-.1 0],[-.1 0],'k-')
ylabel('$\overline{u}_{lin}[m/s]$','interpreter','latex','fontsize',16)
xlabel('$\overline{u}_{meas}[m/s]$','interpreter','latex','fontsize',16)
title('$\overline{u}$','interpreter','latex','fontsize',16)
axis([-.4 .1 -.4 .1])
set(gca,'TickLabelInterpreter','latex')

figure;clf;clear hh
dum = 'osv><^dhp';
for i = 1:length(sve)
  %hhdum = plot(sve(i).meanu,mean(sve(i).F),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on
  % hhdum = plot(mean(sve(i).udum)*abs(mean(sve(i).udum)),mean(sve(i).F/1000),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on
  %hhdum = plot(sve(i).Hrms,mean(sve(i).F/1000),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on
  %plot(sve(i).Hrms,mean(sve(i).F2/1000),['B',dum(i)],'markerfacecolor','k','markersize',6);hold on
  uc = mean(sve(i).udum);
  %hhdum = plot(uc*abs(uc),mean(sve(i).F/1000),['rs'],'markerfacecolor','k','markersize',6);hold on
  %  plot(uc*abs(uc),mean(sve(i).F2/1000),['bs'],'markerfacecolor','k','markersize',6);hold on
  Fd = mean(sve(i).beta*abs(sve(i).udum).*sve(i).udum)*sve(i).hv;
  Feta = sve(i).beta*9.81*(.9*sve(i).Hrms)^3/(6*pi*sve(i).hv);
  hhdum = plot(uc,mean(sve(i).F/1000),['rs'],'markerfacecolor','r','markersize',8);hold on
  hhdum3 = plot(uc,Fd/1000+Feta/1000,['gs'],'markerfacecolor','g','markersize',10);hold on
  %text(uc+.001,mean(sve(i).F/1000),num2str(i))
  hhdum2 =  plot(uc,mean(sve(i).F2/1000),['bs'],'markerfacecolor','b','markersize',8);hold on
  hh(1) = hhdum(1);
  hh(2) = hhdum2;
    hh(3) = hhdum3;
  
  hlabs{i} = dnames(i).name;
end
xl = xlim;

yl = ylim;
plot([xl(1)-.01 .01+xl(2)],0*xl,'k')
plot(0*yl,yl,'k')
hl = legend(hh,'$\beta  d \overline{|u|u} $','$\beta \overline{h |u|u}$','$\beta d \overline{|u|u} + {\bf F}_{\eta}$')
set(hl,'interpreter','latex','fontsize',14,'location','southeast','autoupdate','off')
ylabel('${\bf F}/\rho[m^2/s^2]$','interpreter','latex','fontsize',16)
xlabel('$ \overline{u} [m/s]$','interpreter','latex','fontsize',16)
%xlabel('$u_c/\rho [m^2/s^2]$','interpreter','latex','fontsize',16)
title('Drag vs Current','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter','latex')
%axis([-.4 .1 -.4 .1])
print -dpng dragvscurrent.png

return

dum = 'osv><^dhp';
for i = 1:length(sve)
  figure;clf;clear hh
  %hhdum = plot(sve(i).meanu,mean(sve(i).F),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on
  %hhdum = plot(mean(sve(i).udum)*abs(mean(sve(i).udum)),mean(sve(i).F/1000),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on

  uc = mean(sve(i).udum);
  utilde = sve(i).udum-uc;
  tau = mean(abs(utilde).*utilde);
  t = [0:size(utilde)-1]./100;
  subplot(2,1,1)
  plot(t,utilde,'b-','linewidth',2)
  xlim([0 40])
  xlabel('$t[s]$','interpreter','latex','fontsize',16)
  ylabel('$\tilde{u} [m/s]$','interpreter','latex','fontsize',16)
  title(['',sve(i).name],'interpreter','latex','fontsize',16)
  set(gca,'TickLabelInterpreter','latex')

  subplot(2,1,2)
  histogram(utilde,'Normalization','pdf');hold on
  plot([0 0],[0 1.4],'r','linewidth',3)  
  ylabel('$P$','interpreter','latex','fontsize',16)
  xlabel('$\tilde{u} [m/s]$','interpreter','latex','fontsize',16)
  a = axis;
  text(a(2)-.25*(a(2)-a(1)),a(4)-.22*(a(4)-a(3)),...
       ['$\overline{|\tilde{u}|\tilde{u}}  = ',sprintf('%2.2f',tau),' $'],'interpreter','latex','fontsize',16)

  set(gca,'TickLabelInterpreter','latex')
  %hhdum = plot(uc*abs(uc),mean(sve(i).F/1000),['rs'],'markerfacecolor','k','markersize',6);hold on
  %hh(i) = hhdum(1);
  %hlabs{i} = dnames(i).name;
  dum = get(gcf,'position');
  set(gcf,'position',[dum(1) dum(2) 800 400])
  print('-dpng',['pdf_',sve(i).name,'.png'])
end
%hl = legend(hh,hlabs);
%set(hl,'interpreter','latex','fontsize',14,'location','southeast','autoupdate','off')
%axis([-.4 .1 -.4 .1])



figure;clf;clear hh
for i = 1:length(sve)
  subplot(2,5,i)
  hhdum = plot(std(sve(i).u),sve(i).zu,'ro','markerfacecolor','k','markersize',6);hold on
  hhdum = plot(std(sve(i).w),sve(i).zw,'rs','markerfacecolor','k','markersize',6);hold on
  hhdum = plot(0*std(sve(i).w),sve(i).zw,'k','linewidth',2,'markersize',6);hold on
  axis([-.1 .6 1.4 2])
  if i==1|i==6;ylabel('$z[m]$','interpreter','latex','fontsize',12);end
  if i>5;xlabel('$\sigma_u[m/s]$','interpreter','latex','fontsize',12);end
  set(gca,'TickLabelInterpreter','latex')
end

figure;clf;clear hh
for i = 1:length(sve)
  subplot(2,5,i)
  %utilde = sve(i).u-mean(sve(i).u);
  %wtilde = sve(i).w-mean(sve(i).w);
  P = polyfit(sve(i).zu(1:3),mean(sve(i).utilde(:,1:3).*sve(i).wtilde(:,1:3)),1);
  plot(mean(sve(i).utilde.*sve(i).wtilde),sve(i).zu,'ro','markerfacecolor','k','markersize',6);hold on
  
  plot(polyval(P,sve(i).zu),sve(i).zu,'r-','linewidth',2,'markersize',6);hold on
  plot(0*mean(sve(i).utilde),sve(i).zu,'k-','linewidth',2,'markersize',6);hold on
  plot([-.1 .1],.85+hv*[1 1],'b-','linewidth',2,'markersize',6);hold on
  axis([-.02 .01 .85 2.7])
  if i==1|i==6;ylabel('$z[m]$','interpreter','latex','fontsize',12);end
  if i>5;xlabel('$\overline{\tilde{u} \tilde{w}}$','interpreter','latex','fontsize',12);end
  set(gca,'TickLabelInterpreter','latex')
end


%xlabel('$\overline{u}_{meas}[m/s]$','interpreter','latex','fontsize',16)
%title('$\overline{F},\overline{u}$','interpreter','latex','fontsize',16)



