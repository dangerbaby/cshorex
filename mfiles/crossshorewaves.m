function [wavhyd]=crossshorewaves(in,i,bathy);
  iwaves = 1;
  ibotshear = 1;
  iadvect = 1;
  
  dzbdx = cdiff(in.dx,bathy.zb);  
  xsws = interp_brad(bathy.x,in.swlbc(i)-bathy.zb);
  xmH0 = interp_brad(bathy.x,in.swlbc(i)-.5*in.Hrms(i)-bathy.zb);
  xpH0 = interp_brad(bathy.x,in.swlbc(i)+.5*in.Hrms(i)-bathy.zb);
   tanbeta_foreshore = in.Hrms(i)/(xpH0-xmH0);
   if isnan(xpH0)&~isnan(xsws)
     tanbeta_foreshore = (bathy.zb(end)-in.swlbc(i))/(bathy.x(end)-xsws);
   end    
  % [ib,H0,L0] = surfsim(in.Hrms(i),in.swlbc(i)-bathy.zb(1),in.angle(i),in.Tp(i),tanbeta_foreshore) ;
  %A0= 2.6 + 4.5*ib;
  %A0= 4.5 ;

  if isnan(xsws);xsws=max(in.x);end
  % initialize  
  Hrms = nan(size(bathy.x));Sxx=Hrms;Sxy=Hrms;a=Hrms;Db=Hrms;
  Hm=Hrms;Efr=Hrms;Dr=Hrms;c=Hrms;Er=Hrms;tau_y=Hrms;dtaudv=Hrms;Q=Hrms;k = Hrms;
  eta = max(in.swlbc(i),bathy.zb);    h=eta-bathy.zb;
  [k(1),n(1),c(1)] = dispersion (2*pi/in.Tp(i),h(1));
  %set conditions at offshore boundary
  Hrms(1) = in.Hrms(i);

  Hm(1) = 0.88./k(1).*tanh(in.gamma*k(1).*h(1)/.88);
  Q(1)=probbreaking(Hrms(1),Hm(1));
  Db(1) = dissipation_bore (in.Tp(i),Q(1),Hm(1));
  E(1) = 1/8*9810*Hrms(1)^2;
  Ef(1) = E(1)*c(1)*n(1);
  a(1) = in.angle(i);
  [Sxy(1) Sxx(1) Syy(1)] = radstress (Hrms(1), a(1), n(1), c(1));
  M(1) = in.A0*9.81*h(1)^2;
  Efr(1) = 0;
  Er(1) = Efr(1)/(2*c(1)*cos(a(1)*pi/180));
  Dr(1) = 2*9.81*Er(1)*sin(in.rollerbeta*pi/180)/c(1);
  Ur(1) = undertow_linear(h(1),Hrms(1),in.Tp(i));
  taubx(1)=in.cf(1)*in.Q0^2/h(1)^2;
  [taubx(1) junk] = shear_stress(in.cf(1),Hrms(1),h(1),a(1),in.Tp(i),Ur(1)+in.Q0/h(1),0);
  q2overh(1)=in.Q0^2/h(1); 
  iswash(1)=0;
  %step through transect from offshore to swash
  for j = 2:length(bathy.x)
    if(bathy.x(j)<=xsws) % in the surf zone
      
      %predictor step
      hj = eta(j-1)-bathy.zb(j);
      k(j) = k(j-1);
      n(j) = n(j-1);
      c(j) = c(j-1);
      a(j) = snells(a(1),h(1),in.Tp(i),hj);
      Ef(j) = Ef(j-1)-in.dx*Db(j-1);
      E(j)= max(Ef(j)/(c(j)*n(j)),0);
      Hrms(j) = sqrt(8/9810*E(j));
      [Sxy(j) Sxx(j) Syy(j)] = radstress (Hrms(j), a(j), n(j), c(j));
      %roller
      Efr(j) = Efr(j-1)+in.dx*in.iroll*(Db(j-1)-Dr(j-1));
      Er(j) = Efr(j)/(2*c(j)*cos(a(j)*pi/180));
      Dr(j) = 2*9.81*Er(j)*sin(in.rollerbeta)/c(1);
      Sxy(j)= Sxy(j)+2*Er(j)*cos(a(j)*pi/180)*sin(a(j)*pi/180);
      taubx(j)=in.cf(1)*in.Q0^2/(.5*(hj+h(j-1)))^2;
      q2overh(j)=in.Q0^2/hj; 
      eta(j) = eta(j-1)-2*in.dx/(9810*(h(j-1)+hj))*...
               (iadvect*(q2overh(j)-q2overh(j-1))/in.dx+iwaves*(Sxx(j)-Sxx(j-1))/in.dx+ibotshear*taubx(j));
      eta(j) = max(eta(j),bathy.zb(j));
      h(j) = max(eta(j)-bathy.zb(j),0);
      
      %corrector
      for iter = 1:2
        hj = h(j);
        [k(j),n(j),c(j)] = dispersion (2*pi/in.Tp(i),hj);
        a(j) = snells(a(1),h(1),in.Tp(i),hj);
        Hm(j) = 0.88./k(j).*tanh(in.gamma*k(j).*h(j)/.88);
        Q(j) = probbreaking(Hrms(j),Hm(j));
        Db(j) = dissipation_bore (in.Tp(i),Q(j),Hm(j));
        Ef(j) = Ef(j-1)-in.dx*.5*(Db(j-1)+Db(j));
        E(j)= max(Ef(j)/(c(j)*n(j)),0);
        Hrms(j) = min(Hm(j),sqrt(8/9810*E(j)));
        E(j) = 1/8*9810*Hrms(j)^2; % in case Hm is invoked
        Ef(j) = E(j)*c(j)*n(j);
        Q(j)=probbreaking(Hrms(j),Hm(j));
        %roller
        Efr(j) = Efr(j-1)+in.dx*in.iroll*.5*(Db(j-1)+Db(j)-Dr(j)-Dr(j-1));
        Er(j) = Efr(j)/(2*c(j)*cos(a(j)*pi/180));
        Dr(j) = 2*9.81*Er(j)*sin(in.rollerbeta)/c(j);
        [Sxy(j) Sxx(j) Syy(j)] = radstress (Hrms(j), a(j), n(j), c(j)); % linear representation
                                                                        %Sxx(j)= Sxx(j)+Mr(j)*cos(a(j)*pi/180)^2;
        Sxy(j)= Sxy(j)+2*Er(j)*cos(a(j)*pi/180)*sin(a(j)*pi/180);
        taubx(j)=in.cf(j)*in.Q0^2/(.5*(hj+h(j-1)))^2;
        %[tau_x tau_ypossible(:,j)] = shear_stress(in.cf(j),wavehyd.Hrms(j),wavehyd.h(j),wavehyd.angle(j),in.Tp(i),0,V_possible);
        q2overh(j)=in.Q0^2/hj; 
        eta(j) = eta(j-1) - 2*in.dx/(9810*(h(j-1)+hj)) *...
               (iadvect*(q2overh(j)-q2overh(j-1))/in.dx+iwaves*(Sxx(j)-Sxx(j-1))/in.dx+ibotshear*taubx(j));
        eta(j) = max(eta(j),bathy.zb(j));
        %eta(j) = max(eta(j-1) - 1/(9810*(h(j-1))) * iwaves*(Sxx(j)-Sxx(j-1))/in.dx,zb(j));
        h(j) = eta(j)-bathy.zb(j);
      end
      %M(j) = in.A0*9.81*h(j)^2;
      M(j) = in.A0*9.81*(Hrms(j)/in.gamma)^2;
      iswash(j)=0;
    else
      if iswash(j-1)==0;
        icorrect = 0;
        xintocell = (xsws-bathy.x(j-1))/in.dx;
        dh = icorrect*.45*(bathy.zb(j)-bathy.zb(j-1));% correction factor
        hdeeper = h(j-1)+dh-2*dh*xintocell;
        %hdeeper = min(h(j-1),Hrms(j-1)/in.gamma);        
        %disp(['First swash cell, h(j-1), xintocell = ',num2str(h(j-1)),'  ',num2str(xintocell)])
      else
        hdeeper = h(j-1);
      end
      s = max(0,dzbdx(j));
      loss = (in.cf(j)+(hdeeper>0)*s)*in.dx/(2*in.A0);
      h(j) = max(hdeeper-loss,0);
      Hrms(j) = in.gamma*h(j);
      M(j) = in.A0*9.81*(Hrms(j)/in.gamma)^2;
      eta(j) = bathy.zb(j)+h(j);
      Ef(j) = NaN;
      [k(j),n(j),c(j)] = dispersion (2*pi/in.Tp(i),h(j));
      if h(j)>eps;iswash(j) = 1;else iswash(j)=NaN;end
    end
  end

  % some logical stuff
  Q(iswash==1) = 1;
  Db(isnan(Db)) = 0;
  a(iswash==1) = a(max(find(iswash==0)));
  
  
  %find runup values
  % z_rw = bathy.zb+in.rwh;
  % sig_eta = Hrms/sqrt(8);

  % x1 = interp_brad(bathy.x,eta+sig_eta-z_rw);
  % z1 = interp1(bathy.x,eta+sig_eta,x1);
  % x2 = interp_brad(bathy.x,eta-z_rw);
  % z2 = interp1(bathy.x,eta,x2);
  % x3 = interp_brad (bathy.x,eta-sig_eta-z_rw);
  % z3 = interp1(bathy.x,eta-sig_eta,x3);
  % alternate formulation using max wetted node as runup position
  %[j1]=min(h==0);
  j1=min(find(h==0));
  if ~isempty(j1)
    wavhyd.runup_2p_x = bathy.x(j1);  
    wavhyd.runup_2p  = bathy.zb(j1);  
  else
    wavhyd.runup_2p = NaN;  
    wavhyd.runup_2p_x = NaN;  
  end
  dzbdxdum = 0;
  %wavhyd.runup_mean  = (z1+z2+z3)/3;
  %wavhyd.runup_std   = (z1-z3)/2;
  %wavhyd.runup_13    = wavhyd.runup_mean +(2+0*dzbdxdum)*wavhyd.runup_std;
  %wavhyd.runup_2p    = wavhyd.runup_mean + 1.4*(wavhyd.runup_13-wavhyd.runup_mean);
  %wavhyd.runup_mean_x= interp_brad(bathy.x,wavhyd.runup_mean-bathy.zb);
  %wavhyd.runup_2p_x  = interp_brad(bathy.x,wavhyd.runup_2p-bathy.zb);
  wavhyd.h           = h;
  wavhyd.xsws        = xsws;
  wavhyd.hsws        = interp1(bathy.x,h,xsws);
  wavhyd.Hrms        = Hrms;
  wavhyd.Hm          = Hm;
  wavhyd.Q           = Q;
  wavhyd.Ef          = Ef;
  wavhyd.Efr         = Efr;
  wavhyd.Er          = Er;
  wavhyd.eta         = eta;
  wavhyd.h           = h;
  wavhyd.Db          = Db;
  wavhyd.Dr          = Dr;
  wavhyd.Sxx         = Sxx;
  wavhyd.Sxy         = Sxy;
  wavhyd.angle       = a;
  wavhyd.c           = c;
  wavhyd.k           = k;
  wavhyd.iswash      = iswash;
  wavhyd.A0          = in.A0;
  wavhyd.tanbeta     = tanbeta_foreshore;
  wavhyd.taubx       = taubx;
  

  