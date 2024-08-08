function [in]=set_defaults(in)
for kk = 1:length(in) % loop over in
  if ~isfield(in(kk),'timebc_wave')|isempty(in(kk).timebc_wave);in(kk).timebc_wave = 0;end
  if ~isfield(in(kk),'timebc_surg')|isempty(in(kk).timebc_surg);in(kk).timebc_surg = 0;end
  if ~isfield(in(kk),'iprofl')|isempty(in(kk).iprofl);in(kk).iprofl = 0;end
  if ~isfield(in(kk),'iroll')|isempty(in(kk).iroll);in(kk).iroll = 1;end
  if ~isfield(in(kk),'rollerbeta')|isempty(in(kk).rollerbeta);in(kk).rollerbeta = .1;end
  if ~isfield(in(kk),'cf')|isempty(in(kk).cf);in(kk).cf = in(kk).fw./2;end
  % if ~isfield(in(kk),'xsws')|isempty(in(kk).xsws);
  %   in(kk).xsws = interp_brad(bathy.x,in.swlbc(i)-bathy.zb);
  %   if isnan(in(kk).xsws);in(kk).xsws=max(in(kk).x);end
  % end
  % if ~isfield(in(kk),'tanbeta_foreshore')|isempty(in(kk).tanbeta_foreshore);
  %   xH0 = interp_brad(bathy.x,in.swlbc(i)+in-bathy.zb);
  %   in(kk).tanbeta_foreshore = 5;
  % end
  % if ~isfield(in(kk),'ib')|isempty(in(kk).ib);
  %   [in(kk).ib,H0,L0] = surfsim (in(kk).Hrms,in(kk).swlbc-in(kk).zb(1),in(kk).angle,in(kk).Tp,in(kk).tanbeta_foreshore) ;
  % end
  % if ~isfield(in(kk),'A0')|isempty(in(kk).A0);in(kk).A0= 2.6 + 4.5*in(kk).ib;end
  if ~isfield(in(kk),'A0')|isempty(in(kk).A0);in(kk).A0=3.5;end
  if ~isfield(in(kk),'rwh');in(kk).rwh = .01;end
  if ~isfield(in(kk),'zbhard')|isempty(in(kk).zbhard);in(kk).zbhard = in(kk).zb-100;end
  if ~isfield(in(kk),'verbose');in(kk).verbose = 1;end
  if ~isfield(in(kk),'Q0')|isempty(in(kk).Q0);in(kk).Q0 = 0;end
  if length(in(kk).cf)==1;in(kk).cf=in(kk).cf*ones(size(in(kk).x));end

  
end