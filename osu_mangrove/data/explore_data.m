
clear all
close all
basename = 'HighDensity_h270_hv182_NoWall';
%basename = 'Baseline_h270_hv185_NoWall'
dnames = dir(['./',basename,'/T*']);
%dname = './HighDensity_h270_hv182_NoWall/Trial08/';
%dname = './HighDensity_h270_hv182_NoWall/Trial09/';
cf = .05;
icorrect =1;
cnt = 0;
for j = 1:length(dnames)
  %    for j = 8
  cnt = cnt+1;
  clear p u eta eta_p ubp wbp
  dname = ['./',basename,'/',dnames(j).name,'/'];
  load([dname,'summary.mat'])
  hv = str2num(dname(strfind(dname,'hv')+2:strfind(dname,'hv')+4))/100+.03;
  p = [dat.press.press];
  u = [dat.u.u];u = u(:,2:5);
  xu = [dat.u.x];xu = xu(:,2:5);
  zu = [dat.u.z];zu = zu(:,2:5);
  w = [dat.w.w];w = w(:,2:5);
  zw = [dat.w.z];zw = zw(:,2:5);
  xp = [dat.press.x];

  eta = [dat.wg.eta];
  xwg = [dat.wg.x];
  etaus =[dat.uswg.eta];
  xuswg = [dat.uswg.x];

  activity = mean(std(eta),2);
  startinds = find(mean(abs(eta(:,find(xwg>50))),2)>activity);startind = startinds(10);
  endinds = find(mean(abs(eta(:,find(xwg<50))),2)>activity);endind = endinds(end-10);
  p_init = mean(p(1:round(startind/2),:));
  eta_init = mean(eta(1:round(startind/2),:));
  etaus_init = mean(etaus(1:round(startind/2),:));
  eta = eta(startind:endind,:);
  t = [0:size(eta,1)-1]./100;
  stats = find_stats(t,eta(:,2),4);
  [k,n,c] = dispersion (2*pi/stats.Tp,hv);
  etaus = etaus(startind:endind,:);
  p = p(startind:endind,:);
  u = u(startind:endind,:);
  w = w(startind:endind,:);
  for jj =1:size(u,2)
    ubp(:,jj) = bandpass(u(:,jj),1/100,.5*1/stats.Tp,2*1/stats.Tp);
    wbp(:,jj) = bandpass(w(:,jj),1/100,.5*1/stats.Tp,2*1/stats.Tp);
  end
  
  
  
  utilde = ubp-mean(ubp);
  wtilde = wbp-mean(wbp);
  P = polyfit(zu(1:3),mean(utilde(:,1:3).*wtilde(:,1:3)),1);
  ddzmeanutildewtilde = 1*P(1);
  
  
  for i = 1:size(p,2)
    if dat.press(i).z-dat.press(i).swd<0;
      eta_p(:,i) = p2eta([p(:,i)-mean(p(:,i))],1/100,hv,dat.press(i).z-dat.press(i).swd);
    else
      eta_p(:,i) = dat.press(i).z+p(:,i)/9810;
    end
  end
  udum = u(:,3);
  dx = .1;
  xi = dat.press(1).x:dx:dat.press(6).x;
  Hrmsi = interp1(xwg,sqrt(8)*std(eta),xi);
  ufrometa = eta_p(:,3)*sqrt(9.81/hv);
  ufrometa = ufrometa-mean(ufrometa);
  ufrometa = std(udum)/std(ufrometa)*ufrometa;
  ufrometa = ufrometa + 1*undertow_linear(hv,interp1(xi,Hrmsi,43),stats.Tp);
  taub = 1000*cf*udum.*abs(udum);
  dissb = mean(-taub.*udum);
  De = .041; 
  Cd = 1;
  B = 3.66;L = 18;
  if contains(basename,'Base')
    N = 0;
  elseif contains(basename,'High')
    N= 50*8/(L*B); %number of plants ( and roots) per unit area
  end
  FoverCd = N*1000/2*De*hv*udum.*abs(udum);%force/unit area
  FoverCdp = N*1000/2*De*hv*ufrometa.*abs(ufrometa);%force/unit area
  F2overCd = N*1000/2*De*(hv+eta_p(:,3)).*udum.*abs(udum);% uses p(3) for eta
  F2overCdp = N*1000/2*De*(hv+eta_p(:,3)).*ufrometa.*abs(ufrometa);% uses p(3) for eta
  
  %F2overCd = 1000/2*De*(hv+etaus(:,3)).*u(:,3).*abs(u(:,3));%uses wswg(3) for eta

  %if icorrect;FoverCd=F2overCd;end
  dissvegoverCd = mean(-FoverCd.*udum);
  dissvegoverCdp = mean(-FoverCdp.*ufrometa);
  dissvegoverCd2 = mean(-F2overCd.*udum);
  dissvegoverCd2p = mean(-F2overCdp.*ufrometa);
  Cdexact =((9810*n*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*L)/(dissvegoverCd*L); 
  Cdexactp =((9810*n*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*L)/(dissvegoverCdp*L); 
  Cdexact2 =((9810*n*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*L)/(dissvegoverCd2*L); 
  Cdexact2p =((9810*n*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*L)/(dissvegoverCd2p*L); 
  if N==0;Cdexact = 1;Cdexact2 = 1;end
  %Cd = Cdexact;
  F = N*1000*Cdexact/2*De*hv*udum.*abs(udum);%force per unit area
  Fp = N*1000*Cdexactp/2*De*hv*ufrometa.*abs(ufrometa);%force per unit area
  F2 = N*1000*Cdexact2/2*De*(hv+eta_p(:,3)).*udum.*abs(udum);
  F2p = N*1000*Cdexact2p/2*De*hv*ufrometa.*abs(ufrometa);%force per unit area
  beta = N*1000*Cdexact/2*De;
  beta2 = N*1000*Cdexact2/2*De;
  dissveg = mean(-F.*udum);
  dissvegp = mean(-Fp.*ufrometa);
  dissveg2 = mean(-F2.*udum);
  dissveg2p = mean(-F2p.*ufrometa);


  Efi = (9810/8)*Hrmsi.^2.*n.*c;
  [Sxy Sxx Syy] = radstress (Hrmsi, 0*ones(size(xi)), n*ones(size(xi)), c*ones(size(xi)));
  ddxSxx = cdiff(dx,Sxx);
  %ddxSxx = mean(ddxSxx)*ones(size(ddxSxx));

  etai0 = (mean(p(:,1))-p_init(1))/9810;
  %etai= etai0 - 1/(9810*hv/100)*dx*cumsum(ddxSxx);
  etaiwaves(1) = etai0;
  etaiwavesbotshear(1) = etai0;
  etaiwavesbotshearveg(1) = etai0;
  etaiwavesbotshearvegp(1) = etai0;
  etaiwavesbotshearveg2(1) = etai0;
  etaiwavesbotshearveg2p(1) = etai0;
  etaiwavesbotshearwrs(1) = etai0;%including the wave renyolds stress(wrs)
  etaiwavesbotshearvegwrs(1) = etai0;%including the wave renyolds stress(wrs)
  Efiveg(1)=Efi(1);
  Efiveg2(1)=Efi(1);
  for i = 2:length(xi);

    etaiwaves(i) = etaiwaves(i-1)-(dx/(9810*hv))*ddxSxx(i); 
    etaiwavesbotshear(i)    = etaiwavesbotshear(i-1)-(dx/(9810*hv))*(ddxSxx(i)+mean(taub)); 
    etaiwavesbotshearwrs(i) = etaiwavesbotshearwrs(i-1)-(dx/(9810*hv))*(ddxSxx(i)+1000*hv*ddzmeanutildewtilde+mean(taub)); 
    etaiwavesbotshearveg(i) = etaiwavesbotshearveg(i-1)-(dx/(9810*hv))*(ddxSxx(i)+mean(taub)+mean(F)); 
    etaiwavesbotshearvegwrs(i) = etaiwavesbotshearvegwrs(i-1)-(dx/(9810*hv))*(ddxSxx(i)+1000*hv*ddzmeanutildewtilde+mean(taub)+mean(F)); 
    etaiwavesbotshearvegp(i) = etaiwavesbotshearvegp(i-1)-(dx/(9810*hv))*(ddxSxx(i)+mean(taub)+mean(Fp)); 
    etaiwavesbotshearveg2(i)= etaiwavesbotshearveg2(i-1)-(dx/(9810*hv))*(ddxSxx(i)+mean(taub)+mean(F2)); 
    etaiwavesbotshearveg2p(i) = etaiwavesbotshearveg2p(i-1)-(dx/(9810*hv))*(ddxSxx(i)+mean(taub)+mean(F2p)); 
    Efiveg(i) = Efiveg(i-1)+dx*dissb + dx*dissveg;
    Efiveg2(i) = Efiveg2(i-1)+dx*dissb + dx*dissveg2;
  end
  sve(cnt).u = u;
  sve(cnt).udum = udum;
  sve(cnt).utilde = utilde;
  sve(cnt).wtilde = wtilde;
  sve(cnt).zu = zu;
  sve(cnt).w = w;
  sve(cnt).P = P;
  sve(cnt).zw = zw;
  sve(cnt).ufrometa = ufrometa;
  sve(cnt).meanu = mean(u);
  sve(cnt).ddzmeanutildewtilde = ddzmeanutildewtilde;
  sve(cnt).Tp = stats.Tp;
  [j1 j2]=min(abs(xi-43));
  sve(cnt).Hrms = Hrmsi(j2);
  sve(cnt).hv = hv;
  sve(cnt).undertow = undertow_linear(hv,sve(cnt).Hrms,stats.Tp);
  sve(cnt).F = F;
  sve(cnt).Fp = Fp;
  sve(cnt).F2 = F2;
  sve(cnt).F2p = F2p;
  sve(cnt).stats = stats;
  sve(cnt).name = dnames(j).name;
  sve(cnt).Cd = Cdexact;
  sve(cnt).beta = beta;
  
  figure(1);clf;clear hh hl;
  subplot(2,1,1)
  %plot(xp,sqrt(8)*std(p)/9810,'bs-');hold on
  hh1(1) = plot(xp(1:6),sqrt(8)*std(eta_p(:,1:6)),'bs-','markerfacecolor','b');hold on
  %plot(xuswg,sqrt(8)*std(etaus),'rv','markerfacecolor','r');hold on
  hh1(2) = plot(xwg,sqrt(8)*std(eta),'rs-','markerfacecolor','r');hold on
  hh1(3) = plot(xi,sqrt(8*Efiveg./(9810*n.*c)),'k-','linewidth',2);hold on
  dname2=dname;dname2(dname=='_')='-';
  title([dname2,'  T = ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
    xlabel('$x[m]$','interpreter','latex','fontsize',16)
  ylabel('$H_{rms}[m]$','interpreter','latex','fontsize',16)
  %xlim([10 80])
  a = axis;
  hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2)
  text(37,a(3)+.1*(a(4)-a(3)),'Veg. Section','interpreter','latex','fontsize',14)

  if 0 
    del = .12;
    text(a(1)+.02*(a(2)-a(1)),a(3)+3*del*(a(4)-a(3)),['$C_f ;  C_d  ; C_{d2} = $',...
                        sprintf('%4.2f   %4.2f  %4.2f  ',[cf Cdexact Cdexact2])],'interpreter','latex','fontsize',12)
    %  text(a(1)+.02*(a(2)-a(1)),a(3)+3*del*(a(4)-a(3)),['$\rho \sigma_u^3 ',...
    %                    sprintf('%2.2f ',[1000*std(u(:,3))^3]),'$'],'interpreter','latex','fontsize',12)
    text(a(1)+.02*(a(2)-a(1)),a(3)+2*del*(a(4)-a(3)),['$D_{bot} ; D_{veg} = $',...
                        sprintf('%2.2f  %2.2f',[dissb,dissveg])],'interpreter','latex','fontsize',12)
    text(a(1)+.02*(a(2)-a(1)),a(3)+del*(a(4)-a(3)),['${\bf F} \; {\bf F}(z) = $',...
                        sprintf('%2.2f  %2.2f ',[mean(F) mean(F2)])],'interpreter','latex','fontsize',12)
  else
    
    hl1 = legend(hh1,'P-Guage','C-Guage','Calibrated Model','AutoUpdate','off');  
    set(hl1,'interpreter','latex','fontsize',12,'location','southwest')
   set(gca,'xticklabels',[])
  end
  
  set(gca,'TickLabelInterpreter','latex')

  subplot(2,1,2);
  plot(xp(1:6),1000*(mean(p(:,1:6))-p_init(1:6))/9810,'bs-','markerfacecolor','b');hold on
  plot(xwg,1000*(mean(eta)-eta_init),'rs-','markerfacecolor','r');hold on
  %plot(xuswg,1000*(mean(etaus)-etaus_init),'rv','markerfacecolor','r');hold on
  %plot(xi,1000*etaiwaves,'k-','linewidth',1)
  hh(1) = plot(xi,1000*etaiwavesbotshear,'k--','linewidth',1);
  %hh(2) = plot(xi,1000*etaiwavesbotshearveg,'k-.','linewidth',2);
  %plot(xi,1000*etaiwavesbotshearwrs,'m--','linewidth',2);
  %plot(xi,1000*etaiwavesbotshearveg2p,'m-','linewidth',2);
  hh(2) = plot(xi,1000*etaiwavesbotshearveg2,'k-','linewidth',2);
  hl = legend(hh,'$S_{xx}, \tau_b$','$S_{xx}, \tau_b$, {\bf F}','$S_{xx}, \tau_b, {\bf F}(z)$');
  hl = legend(hh,'$S_{xx}, \tau_b$','$S_{xx}, \tau_b, {\bf F}$','AutoUpdate','off');
  set(hl,'interpreter','latex','fontsize',14,'location','northwest')
  dum = ylim;
  if dum(2)-dum(1)<100
    ylim([mean(ylim)-40 mean(ylim)+40])
  end
  xlabel('$x[m]$','interpreter','latex','fontsize',16)
  ylabel('$\overline{\eta}[mm]$','interpreter','latex','fontsize',14)
  set(gca,'TickLabelInterpreter','latex')
  %print('-dpng',[Data(ind).name(1:end-4),'.png'])
  dname2=dname;dname2(dname=='_'|dname=='/')='-';
  a = axis;
  hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2)
  % text(a(1)+.75*(a(2)-a(1)),a(3)+.3*(a(4)-a(3)),...
  %      ['$d \overline{|u| u} = ',sprintf('%2.2f',mean(F/beta)),'m^3/s^2$'],'interpreter','latex','fontsize',12)
  % text(a(1)+.75*(a(2)-a(1)),a(3)+.1*(a(4)-a(3)),...
  %      ['$\overline{h|u| u} - d \overline{|u| u}$ = ',sprintf('%2.2f',mean(F2/beta2) -mean(F/beta))],'interpreter','latex','fontsize',12)
  % dumstr = {['$\overline{u} = ',sprintf('%2.2f',mean(udum)),'\frac{m}{s}$'];
  %           ['$d \overline{|u| u} = ',sprintf('%2.2f',mean(F/beta)),'\frac{m^3}{s^2}$'];
  %           ['$\overline{h|u| u} - d \overline{|u| u} = ',sprintf('%2.2f',mean(F2/beta2) -mean(F/beta)),'\frac{m^3}{s^2}$']}
  dt = .01;
  tdum = 0:dt:100*sve(cnt).Tp;
  umean = mean(udum);
  utilde = udum-mean(udum);
  usym = umean + sqrt(2)*std(udum)*sin(tdum*2*pi/sve(cnt).Tp);
  Fc = beta*hv*mean(abs(usym).*usym);
  Fasym =.5*sve(cnt).hv*sve(cnt).beta*(mean(abs(udum).*udum)-mean(abs(umean-utilde).*(umean-utilde)));
  
  
  dumstr = {['$F_c = ',sprintf('%2.1f',Fc),' \frac{kg}{m s^2}$'];
            ['$F_{asym} = ',sprintf('%2.1f',mean(Fasym)),' \frac{kg}{m s^2}$'];
            ['$F_{FS} = ',sprintf('%2.1f',mean(F2) -mean(F)),' \frac{kg}{m s^2}$']}
  ha = annotation('textbox', [0.71, 0.14, 0.19, 0.14], 'String',dumstr,'interpreter','latex',...
                  'fontsize',10,'backgroundcolor','w')
  print('-dpng','-r300',[dname2(3:end-1),'.png'])
end


figure(2);clf;clear hh
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

figure(3);clf;clear hh
dum = 'osv><^dhp';
for i = 1:length(sve)
  %hhdum = plot(sve(i).meanu,mean(sve(i).F),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on
  hhdum = plot(1000*mean(sve(i).udum)*abs(mean(sve(i).udum)),mean(sve(i).F),['r',dum(i)],'markerfacecolor','k','markersize',6);hold on
  hh(i) = hhdum(1);
  hlabs{i} = dnames(i).name;
end
%hl = legend(hh,hlabs);
xl = xlim;
plot(xl,0*xl,'k')
set(hl,'interpreter','latex','fontsize',14,'location','northwest','autoupdate','off')
ylabel('$\overline{F}[\rho m^2/s^2]$','interpreter','latex','fontsize',16)
xlabel('$\rho \overline{|u|}\; \overline{u} [m^2/s^2]$','interpreter','latex','fontsize',16)
title('Drag vs Current','interpreter','latex','fontsize',16)
set(gca,'TickLabelInterpreter','latex')
%axis([-.4 .1 -.4 .1])


figure(4);clf;clear hh
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

figure(5);clf;clear hh
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



