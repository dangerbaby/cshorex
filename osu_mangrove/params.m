

  in.name   = name;
  in.iprofl = 0;          % 0 = no morph, 1 = run morph, 1.1 = run morph without initial smoothing
  in.isedav = 0;          % 0 = unlimited sand, 1 = hard bottom
  in.iroll  = 0;          % 0 = no roller, 1 = roller
  in.iwind  = 0;          % 0 = no wind effect
  in.dx     = 1;          % constant dx 
  in.gamma  = .6;         % shallow water ratio of wave height to water depth
  in.sporo  = 0.4;        % sediment porosity                        
  in.d50 = .25;            % d_50 in mm
  in.wf = vfall(in.d50,20,0); % fall velocity
  in.sg = 2.65;           % specific gravity
  in.effb   = 0.002;      % suspension efficiency due to breaking eB     
  in.efff   = 0.004;       % suspension efficiency due to friction ef 
  in.slp    = .5;         % suspended load parameter               
  in.slpot  = .1;         % overtopping suspended load parameter               
  in.tanphi = .630;       % tangent (sediment friction angle)        
  in.blp    = 0.001;      % bedload parameter                        
                          %in.rwh = .02;           % numerical rununp wire height 
  in.cf = .03;     % bottom friction factor
                     %in.A0 = 5;
  in.nut = 0.;

  % boundary conditions and timing
  dt = 1800;         % time interval in seconds for wave and water level conditions
  ftime = .1*dt;      % [sec] final time, dictates model duration
  in.timebc_wave = [0:dt:ftime];

  in.timebc_surg = in.timebc_wave;
  %in.nwave =  length(in.timebc_wave); in.nsurg = in.nwave;
  dum = ones(1,length(in.timebc_wave));
  in.Tp= 8*dum;        % constant spectral peak period in seconds
  in.Hrms = .8+0*cumsum(dum);
  in.Wsetup = 0*dum;   % wave setup at seaward boundary in meters
  in.swlbc = -0+.00*cumsum(dum);% water level at seaward boundary in meters
  in.angle = 0*10*dum;    % constant incident wave angle at seaward boundary in

  % Idealized numerical tank
  Lx = 200;              % length of domain
  zb_off = -3;           % offshore bottom position (should be negative)
  zb_on = 1.5;             % onshore bottom position (should be pos)
  flat_length = 0;     % length of flat portion at seaward end of numerical tank
  x = [0 flat_length Lx];% x points
  zb = [zb_off zb_off zb_on]; % zb points
  in.x = 0:in.dx:Lx;
  [j1 j2] = unique(x); 
  in.zb = interp1(x(j2),zb(j2),in.x);
  dunecrestx = .6*Lx;
  dunestd = sqrt(.05*Lx);
  dunetoexland = dunecrestx+1*dunestd;

  in.zb = in.zb + 0*.3*exp(-(in.x-dunecrestx).^2/(dunestd)^2);
  dunetoeland= interp1(in.x,in.zb,dunetoexland);
  
  %in.zb(in.x>dunetoexland) = dunetoeland;
   
  
  %in.zbhard = in.zb;in.zbhard(find(in.x<100))= in.zbhard(find(in.x<100))-1;
  %in.zb = window(in.zb,11);
  %in.fw = in.fric_fac*ones(size(in.zb)); % cross-shore values of bot fric
  in.cf = in.cf*ones(size(in.zb)); % cross-shore values of bot fric

  in.Q0 = 0;


