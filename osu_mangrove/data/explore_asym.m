dt = .01;


for i = 1:5
  tdum = 0:dt:100*sve(i).Tp;
  
  u = sve(i).udum;
  umean = mean(u);
  utilde = u-mean(u);
  usym = umean + sqrt(2)*std(u)*sin(t*2*pi/sve(i).Tp);
  
  
  
  % disp(num2str([abs(umean)*umean mean(abs(u).*u) mean(abs(utilde).*utilde)]))
  disp([sve(i).name,' F_current = ',num2str(sve(i).hv*sve(i).beta*mean(abs(usym).*usym))])
  disp([sve(i).name,' F_asym = ',num2str(.5*sve(i).hv*sve(i).beta*(mean(abs(u).*u)-mean(abs(umean-utilde).*(umean-utilde))))])

  %disp(num2str([mean(abs(u).*u) mean(abs(umean-utilde).*(umean-utilde))  mean(abs(usym).*usym) ]))
 
end
 
 