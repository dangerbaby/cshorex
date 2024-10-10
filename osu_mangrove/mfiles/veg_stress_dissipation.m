%function [out] = veg_stress_dissipation(in)
%in.cd         veg drag coeff 
%in.vegheight [m] vegitation height from bottom 
%in.Hrms [m]     wave height
%in.d [m]       meanwater depth
%in.alpha [deg] wave angle rel. to x
%in.T [sec]       wave period
%in.U [m/s]       Current in x
%in.V [m/s]       Current in y  
%optional:
%in.inonlin = 0 for linear waves, 1 for nonlin waves
%out.eta[m]  time-depdendent free-surface position
%out.u  [m/s] total near-bed vel in x
%out.v  [m/s] total near-bed vel in y
%out.fx [rho m/s^2)] time, depth dependent force in x
%out.Fx [rho m^2/s^2)] time-dependent depth-integrated force in x
%out.meanFx [rho m^2/s^2)] phase-averged depth-integrated force in x
%note rho = 1000 is used herein
%
function [out] = veg_stress_dissipation(in)
if ~isfield(in,'V');in.V=0;end
if ~isfield(in,'alpha');in.alpha=0;end
if ~isfield(in,'vegheight');in.vegheight=Inf;end
numz = 21;
numt = 51;
rho = 1000;
out.t = linspace(0,in.T,numt);
w = 2*pi/in.T;
[k,n,c] = dispersion (w,in.d);
out.z = repmat(linspace(0,1,numz)',1,numt);

if ~isfield(in,'inonlin')|in.inonlin==0
  out.eta  = (in.Hrms./2)*sin(w*out.t);
else
  
  in2.t = out.t;
  in2.T = in.T;
  in2.Hs = sqrt(2)*in.Hrms;
  in2.k  = k;
  in2.d = in.d;
  [in2.Ur] = 1*ursell(in2.Hs,in2.k,in2.d);
  [out2] = nonlin_free_surface_shape (in2);
  out.eta = out2.eta_t;
  %wave_vel = (2*pi/T)*out2.eta_t./sinh(k*h);
end

out.h = in.d+out.eta;
out.z = out.z.*(in.d+repmat(out.eta,numz,1));% dist from bed
dz = out.z(2,:)-out.z(1,:);
out.wave_vel = w*out.eta.*cosh(k*out.z)./sinh(k*in.d);
out.u = in.U+out.wave_vel*cos(in.alpha*pi/180);
out.v = in.V+out.wave_vel*sin(in.alpha*pi/180);
abs_vel = sqrt(out.u.^2+out.v.^2);
out.fx = 0.5*rho*in.cd*in.N*in.D*abs_vel.*out.u;
%out.fx = out.fx.*(out.z<=in.vegheight);
out.Fx = sum(out.fx).*dz;
out.meanFx = mean(out.Fx);
dx = -0.5*rho*in.cd*in.N*in.D*abs_vel.^2.*abs(out.u);
out.dx = dx.*(out.z<=in.vegheight);
out.Dx = sum(out.dx).*dz;
out.meanDx = mean(out.Dx);




return

u = U+wave_vel*cos(alpha*pi/180);
v = V+wave_vel*sin(alpha*pi/180);
abs_vel = sqrt(u.^2+v.^2);
out.tau_x = rho*cf*(abs_vel.*u);
out.tau_y = rho*cf*(abs_vel.*v);
out.meantau_x = mean(out.tau_x); 
out.meantau_y = mean(out.tau_y); 
out.eta = eta;
out.t = t;
out.u = u;
out.v = v;


