[j1 j2] = min(abs(xi-43));


figure(1);clf
plot(t,udum,'r');hold all
plot(sav(j2).t+3.8,sav(j2).u,'b')
ylabel('$u $','interpreter','latex','fontsize',16)
xlim([0 20]);



figure(2);clf
plot(t,hv+eta_p(:,3),'r');hold all
plot(sav(j2).t+3.8,sav(j2).h,'b')
ylabel('$h $','interpreter','latex','fontsize',16)
xlim([0 20]);

figure(3);clf
plot(t,udum.*abs(udum).*(hv+eta_p(:,3)),'r');hold all
plot(sav(j2).t+3.8,mean(sav(j2).u).*abs(mean(sav(j2).u)).*sav(j2).h,'b','linewidth',2)
ylabel('$h u |u| $','interpreter','latex','fontsize',16)
xlim([0 20]);



figure(4);clf
plot(t,0*F2,'k');hold all
hf(1) = plot(t,F2,'r');hold all
plot(t,mean(F2)+0*t,'r');hold all
hf(2) = plot(sav(j2).t+3.8,sav(j2).Fx,'b','linewidth',2)
plot(sav(j2).t+3.8,sav(j2).meanFx+0*sav(j2).t,'b','linewidth',2)
ylabel('$F_x $','interpreter','latex','fontsize',16)
xlim([0 20]);
title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
hl = legend(hf,'Using Measured Data','Using Modeled Data','AutoUpdate','off');
set(hl,'interpreter','latex','fontsize',14,'location','northwest')
