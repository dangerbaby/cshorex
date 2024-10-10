clear all
close all
addpath('./mfiles');
basename = 'HighDensity_h270_hv182_NoWall';
%basename = 'Baseline_h270_hv185_NoWall'
%dnames = dir(['~/data/osu_mangrove/',basename,'/T*']);
dnames = dir(['./data/',basename,'/T*']);
%dname = './HighDensity_h270_hv182_NoWall/Trial08/';
%dname = './HighDensity_h270_hv182_NoWall/Trial09/';
cf = .05;
icorrect =1;
cnt = 0;
offset = [3. 3.9 2.8 3.5 3.8 3.8]
%for j = 1:length(dnames)
for j = 1:6
  cnt = cnt+1;
  clear p u eta eta_p ubp wbp
  dname = ['~/data/osu_mangrove/',basename,'/',dnames(j).name,'/'];
  dname = ['./data/',basename,'/',dnames(j).name,'/'];
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

  for i = 1:size(p,2)
    if dat.press(i).z-dat.press(i).swd<0;
      eta_p(:,i) = p2eta([p(:,i)-mean(p(:,i))],1/100,hv,dat.press(i).z-dat.press(i).swd);
    else
      eta_p(:,i) = dat.press(i).z+p(:,i)/9810;
    end
  end

  % get exact value of Cd
  dx = .1;
  xi = dat.press(1).x:dx:dat.press(6).x;
  Hrmsi = interp1(xwg,sqrt(8)*std(eta),xi);
  %Hrmsi = interp1(xp,sqrt(8)*std(eta_p),xi);
  [Sxy Sxx Syy] = radstress (Hrmsi, 0*ones(size(xi)), n*ones(size(xi)), c*ones(size(xi)));
  ddxSxx = cdiff(dx,Sxx);
  udum = u(:,3);
  %udum = udum-mean(udum);
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
  F2overCd = N*1000/2*De*(hv+eta_p(:,3)).*udum.*abs(udum);% uses p(3) for eta
  dissvegoverCd2 = mean(-F2overCd.*udum);
  Cdexact2 =((9810*n*c/8)*(Hrmsi(end)^2-Hrmsi(1)^2) - dissb*L)/(dissvegoverCd2*L); 
  if N==0;Cdexact = 1;Cdexact2 = 1;end
  %Cd = Cdexact;
  F2 = N*1000*Cdexact2/2*De*(hv+eta_p(:,3)).*udum.*abs(udum);
  d = -(abs(F2).*abs(udum));
  dissveg2 = -mean(abs(F2).*abs(udum));
  
  modelEf(1) = (9810*n*c/8)*Hrmsi(1)^2;
  modelHrms(1) = sqrt(8*modelEf(1)/(9810*n*c));
  datEf(1) = (9810*n*c/8)*Hrmsi(1)^2;
  datHrms(1) = sqrt(8*modelEf(1)/(9810*n*c));
  in.cd = 1*Cdexact2;    in.d = hv;
  in.T = stats.Tp;  in.N = N;  in.D = De;  in.alpha = 0;  in.vegheight = Inf;
  %in.U = undertow_linear(in.d,in.Hrms,in.T);
  in.U = mean(udum);
  in.Hrms = 1.*modelHrms(1);
  in.inonlin = 1;
  [out] = veg_stress_dissipation(in);
  sav(1) = out;
  modeldissveg(1) = out.meanDx;
  modeleta(1) = (mean(p(:,1))-p_init(1))/9810;
  modelh(1) = hv+modeleta(1);
  dateta(1) = (mean(p(:,1))-p_init(1))/9810;
  eta0a(1) = (mean(p(:,1))-p_init(1))/9810;
  eta0b(1) = (mean(p(:,1))-p_init(1))/9810;
  [Sxy Sxx(1) Syy] = radstress(modelHrms(1),0, n, c);
  for k = 2:length(xi)
    datEf(k) = datEf(k-1)+dx*(dissb+dissveg2);
    modelEf(k) = modelEf(k-1)+dx*(dissb+modeldissveg(k-1));
    modelHrms(k) = sqrt(8*modelEf(k)/(9810*n*c));
    datHrms(k) = sqrt(8*datEf(k)/(9810*n*c));
    in.Hrms = modelHrms(k);
    [out] = veg_stress_dissipation(in);
    sav(k) = out;
    modeldissveg(k) = out.meanDx;
    modelF(k) = out.meanFx;
    [Sxy Sxx(k) Syy] = radstress(modelHrms(k),0, n, c);
    modeleta(k)= modeleta(k-1)-1*(Sxx(k)-Sxx(k-1))/(9810*hv)  - (dx/(9810*hv))*(mean(taub)+modelF(k)); 
    dateta(k)= dateta(k-1)  - (dx/(9810*hv))*(ddxSxx(k)+mean(taub)+mean(F2)); 
    eta0a(k)= eta0a(k-1)-1*(Sxx(k)-Sxx(k-1))/(9810*hv)  - (dx/(9810*hv))*(mean(taub)+0*modelF(k)); 
    eta0b(k)= eta0b(k-1)  - (dx/(9810*hv))*(ddxSxx(k)+mean(taub)+0*mean(F2)); 
  end
  
  
  
  
  
  figure(cnt);clf;clear hh hl;
  subplot(2,1,1)
  %plot(xp,sqrt(8)*std(p)/9810,'bs-');hold on
  hh1(1) = plot(xp(1:6),sqrt(8)*std(eta_p(:,1:6)),'bs-','markerfacecolor','b');hold on
  %plot(xuswg,sqrt(8)*std(etaus),'rv','markerfacecolor','r');hold on
  hh1(2) = plot(xwg,sqrt(8)*std(eta),'rs-','markerfacecolor','r');hold on
  hh1(3) = plot(xi,modelHrms,'k-','linewidth',2);hold on
  hh1(3) = plot(xi,datHrms,'k--','linewidth',2);hold on
  %dname2=dname;dname2(dname=='_')='-';
  title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  %xlabel('$x[m]$','interpreter','latex','fontsize',16)
  ylabel('$H_{rms}[m]$','interpreter','latex','fontsize',16)
  a = axis;
  axis manual
  hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2)
  text(37,a(3)+.1*(a(4)-a(3)),'Veg. Section','interpreter','latex','fontsize',14)
  text(15,a(3)+.1*(a(4)-a(3)),['$C_d = ',sprintf('%2.2f',Cdexact2),'$'],'interpreter','latex','fontsize',14)

  subplot(2,1,2);
  plot(xp(1:6),1000*(mean(p(:,1:6))-p_init(1:6))/9810,'bs-','markerfacecolor','b');hold on
  plot(xwg,1000*(mean(eta)-eta_init),'rs-','markerfacecolor','r');hold on
  %plot(xuswg,1000*(mean(etaus)-etaus_init),'rv','markerfacecolor','r');hold on
  hh(1) = plot(xi,1000*modeleta,'k-','linewidth',2);
  hh(2) = plot(xi,1000*dateta,'k--','linewidth',2);
  hh(3) = plot(xi,1000*eta0a,'g-','linewidth',2);
  %  hh(2) = plot(xi,1000*eta0b,'g--','linewidth',2);
  %hh(2) = plot(xi,1000*etaiwavesbotshearveg,'k-.','linewidth',2);
  %plot(xi,1000*etaiwavesbotshearwrs,'m--','linewidth',2);
  %plot(xi,1000*etaiwavesbotshearveg2p,'m-','linewidth',2);
  %hh(2) = plot(xi,1000*etaiwavesbotshearveg2,'k-','linewidth',3);
  hl = legend(hh,'Modeled $ \eta, u $','Measued $\eta, u $','${\bf F}= 0$','AutoUpdate','off');
  %hl = legend(hh,'$S_{xx}, \tau_b$','$S_{xx}, \tau_b, {\bf F}$','AutoUpdate','off');
  set(hl,'interpreter','latex','fontsize',14,'location','northwest')
  dumyl = ylim;
  if dumyl(2)-dumyl(1)<100
    ylim([mean(ylim)-40 mean(ylim)+40])
  end
  xlabel('$x[m]$','interpreter','latex','fontsize',16)
  ylabel('$\overline{\eta}[mm]$','interpreter','latex','fontsize',14)
  set(gca,'TickLabelInterpreter','latex')

  dname2=dname;dname2(dname=='_'|dname=='/')='-';
  print('-dpng',[dname2(21:end-1),'.png'])
  a = axis;
  hf=fill([36 54 54 36],[a(3) a(3) a(4) a(4)],[.2 .5 .2]);set(hf,'facealpha',.2)
  axis(a)
  
  [j1 j2] = min(abs(xi-43));
  figure(2);clf
  subplot(3,1,1)
  plot(t,udum,'r','linewidth',2);hold all
    plot(t,0*udum,'k','linewidth',1);hold all
  plot(sav(j2).t+offset(j),mean(sav(j2).u),'b','linewidth',2)
  ylabel('$u$','interpreter','latex','fontsize',16)
  %xlabel('$t[s]$','interpreter','latex','fontsize',16)
  xlim([0 20]);
  title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  %print('-dpng',[dname2(21:end-1),'u.png'])

  subplot(3,1,2)
  plot(t,-d,'r','linewidth',2);hold all
  plot(t,0*t+mean(-d),'r--','linewidth',1);hold all
  plot(sav(j2).t+offset(j),-sav(j2).Dx,'b','linewidth',2)
  plot(sav(j2).t+offset(j),ones(size(sav(j2).Dx))*mean(-sav(j2).Dx),'b--','linewidth',1)
  ylabel('$-D$','interpreter','latex','fontsize',16)
  %xlabel('$t[s]$','interpreter','latex','fontsize',16)
  xlim([0 20]);
  %title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  %print('-dpng',[dname2(21:end-1),'u.png'])
  
  subplot(3,1,3)
  plot(t,0*F2,'k');hold all
  hf(1) = plot(t,F2,'r','linewidth',2);hold all
  plot(t,mean(F2)+0*t,'r--');hold all
  hf(2) = plot(sav(j2).t+offset(j),sav(j2).Fx,'b','linewidth',2)
  plot(sav(j2).t+offset(j),sav(j2).meanFx+0*sav(j2).t,'b--','linewidth',1)
  ylabel('$F_x $','interpreter','latex','fontsize',16)
  xlabel('$t[s]$','interpreter','latex','fontsize',16)
  xlim([0 20]);
  %title([strrep(dname(21:end),'_','-'),' ',sprintf('%2.2f',Hrmsi(1)),' ',sprintf('%2.2f',stats.Tp)],'interpreter','latex','fontsize',12)
  hl = legend(hf,'Using Measured Data','Using Modeled Data','AutoUpdate','off');
  set(hl,'interpreter','latex','fontsize',14,'location','northeast')
  print('-dpng',[dname2(21:end-1),'_u_d_F.png'])
end

