function [out] = nonlin_free_surface_shape (in)
%function [out] = nonlin_free_surface_shape (in)
% required
% in.T  = Period[s]
% in.Ur  = Ursell num
% in.Hs  = Sig wave height [m]
% in.k   = wavnumber [1/m]
% in.d   = water depth [m]
% optional:
% in.x   = cross-shore position (use closely-spaced vector ) [m]
% in.t   = time series for out
% out.eta_x = spatial variation of eta in independet columns;
% out.eta_t = time variation of eta in independet columns;
% out.eta_xt = time, space variation of eta with time as columns;

in.Ur=reshape(in.Ur,1,[]);
in.Hs=reshape(in.Hs,1,[]);
in.k=reshape(in.k,1,[]);
in.d=reshape(in.d,1,[]);

p=[0 0.857 -0.417 0.297 0.815 0.672];
B = p(1)+(p(2)-p(1))./(1+exp((p(3)-log10(in.Ur))/p(4)));
psi = -90 +90*tanh(p(5)./in.Ur.^p(6));
out.Sk = B.*cos(psi*pi/180);
out.As = B.*sin(psi*pi/180);
rdum = [0:.01:.99]';
b = rdum./(1+sqrt(1-rdum.^2));
Bdum = 3*b./(2*(1-b.^2));
phi = -psi - 90;
r = interp1(Bdum,rdum,B);
if ~isfield(in,'t')
  t = linspace(0,in.T,101)';t = repmat(t(1:end-1),1,length(in.Ur));numt = size(t,1);
else
  t = in.t;numt = size(t,1);
end


if ~isfield(in,'x')
  x = linspace(0,1,91)';x = repmat(x(1:end-1),1,length(in.Ur)).*repmat(2*pi./in.k,length(x)-1,1);numx = size(x,1);

  % function  of time with independent columns:
  kx = 0;
  omegat = t*2*pi/in.T;
  Hs2 = repmat(in.Hs,numt,1);
  r2 = repmat(r,numt,1);
  phi2 = repmat(phi,numt,1)*pi/180;
  out.eta_t = (sin(-kx+omegat)+r2.*sin(phi2)./(1+sqrt(1-r2.^2)))./(1-r2.*cos(-kx+omegat+phi2));
  out.eta_t = Hs2./(4*repmat(std(out.eta_t),numt,1)).*out.eta_t;
  out.t = t;

  % function of x with indepenednt columns:
  kx = repmat(in.k,numx,1).*x ;
  omegat = 0;
  Hs2 = repmat(in.Hs,numx,1);
  r2 = repmat(r,numx,1);
  phi2 = repmat(phi,numx,1)*pi/180;
  out.eta_x = (sin(-kx+omegat)+r2.*sin(phi2)./(1+sqrt(1-r2.^2)))./(1-r2.*cos(-kx+omegat+phi2));
  out.eta_x = Hs2./(4*repmat(std(out.eta_x),numx,1)).*out.eta_x;
  out.x = x;

else
  dx = in.x(2)-in.x(1);
  in.x=reshape(in.x,1,[]);
  
  out.x = repmat(in.x,numt,1);
  out.t = t;
  kx = repmat(dx*cumsum([0 in.k(1:end-1)]),numt,1);
  omegat = t*2*pi/in.T;
  Hs2 = repmat(in.Hs,numt,1);
  r2 = repmat(r,numt,1);
  phi2 = repmat(phi,numt,1)*pi/180;
  out.eta_xt = (sin(-kx+omegat)+r2.*sin(phi2)./(1+sqrt(1-r2.^2)))./(1-r2.*cos(-kx+omegat+phi2));
  out.eta_xt = Hs2./(4*repmat(std(out.eta_xt),numt,1)).*out.eta_xt;
  
  
end
