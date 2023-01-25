function [in]=set_defaults(in)
for kk = 1:length(in) % loop over in

  if ~isfield(in(kk),'iprofl')|isempty(in(kk).iprofl);in(kk).iprofl = 0;end
  if ~isfield(in(kk),'iroll')|isempty(in(kk).iroll);in(kk).iroll = 1;end
  if ~isfield(in(kk),'rollerbeta')|isempty(in(kk).rollerbeta);in(kk).rollerbeta = .1;end
  if ~isfield(in(kk),'cf');in(kk).cf = in(kk).fw./2;end
  if ~isfield(in(kk),'A0')|isempty(in(kk).A0);in(kk).A0 = 5;end
  if ~isfield(in(kk),'rwh');in(kk).rwh = .01;end
  if ~isfield(in(kk),'zbhard')|isempty(in(kk).zbhard);in(kk).zbhard = in(kk).zb-100;end
  if ~isfield(in(kk),'verbose');in(kk).verbose = 1;end
  if length(in(kk).cf)==1;in(kk).cf=in(kk).cf*ones(size(in(kk).x));end
  
  
end