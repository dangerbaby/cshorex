% function [x2] = bandpass(x,dt,fcutlow,fcuthigh)
% 
function [x2] = bandpass(x,dt,fcutlow,fcuthigh)
N = size(x,1);
df = 1/(N*dt);
f=[ 0:N/2  N/2-1:-1:1]'.*df;
C = fft(x);
ind = find(f>fcuthigh);C(ind,:)=0;
ind = find(f<fcutlow);C(ind,:)=0;
x2 = real(ifft(C));
