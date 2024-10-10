% [eta] = p2eta(p,dt,d,z,f_max)
% p[kg/(m s^2)] = dynamic pressure time series with mean(p) = 0
% dt[s] = 1/sample frequency
% d[m] = still-water or slowly varied depth
% z[m] = position of the pressure sensor relative to the 
%        mean water level.  z must be negative.
% f_max = maximum frequency to include in the p->eta conversion
%
function [eta] = p2eta(p,dt,d,z,f_max)
if d<=0;
  error('Error: d must be finite positive ')
end
if z>=0;
  error('Error: z must be negative to be located below the MWL ')
end
if mean(p)>=std(p)/100;
  disp(['Warning: mean(p) = ',num2str(mean(p))])
  disp(['std(p) = ',num2str(std(p))])
  disp('A pressure record with the mean removed will be used in this analysis')
  p=p-mean(p);
  %  error('Error: p is dynamic pressure, and mean(p) must equal zero')
end
p=p(:);
flag=0;
if mod(length(p),2)==0;
  disp('Warning: p2eta expects an odd number of data points.')
  disp('The data is padded one point for this analysis')
  p=[p; p(end)];
  flag=1;
end

N = length(p);
fny = 1/(2*dt);
T = dt*(length(p)-1);omega1 = 2*pi/T;f1 = 1/T;
Cn=fft(p);
Cn2=2*Cn(2:(length(Cn)-1)/2+1);
f = [f1:f1:fny]';
k= dispersion(2*pi*f,d);
kp = cosh(k*(d+z))./cosh(k*d);
%kp = ones(size(kp));

if exist('f_max')
  ind_bad = find(f>f_max);
else
  ind_bad = find(kp<.1|isnan(kp));
end
disp(['Note: all frequencies greater than ',num2str(min(f(ind_bad))),'[1/s]'])
disp(['     (T<',num2str(1/min(f(ind_bad))),') are neglected in this analysis'])
%whos
kp2 = [0 ;kp; flipud(kp)];
f2 = [Inf ;f; flipud(f)];

if exist('f_max')
  ind_bad = find(f2>f_max);
else
  ind_bad = find(kp2<.2|isnan(kp2));
end
kp2(ind_bad)=1;
tot = (sum(abs(Cn/N).^2)/2);loss =(sum(abs(Cn(ind_bad)/N).^2)/2);
Cntemp = (Cn./kp2);
Cntemp(ind_bad)=0;
%disp([num2str(100*loss/tot),'% of the variance in pressure is not represented in eta']) 
eta= real(ifft(Cntemp/9810));
if flag
  eta = eta(1:end-1);
end
