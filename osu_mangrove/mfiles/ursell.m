%function [Ur] = ursell (Hmo,k,h)
function [Ur] = ursell (Hmo,k,h)
Hmo=reshape(Hmo,1,[]);
k=reshape(k,1,[]);
h=reshape(h,1,[]);
Ur = (3/8)*(Hmo.*k)./((k.*h).^3);
