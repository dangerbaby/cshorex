function check_in(in)
for kk = 1:length(in) % loop over in
  fns =  fieldnames(in(kk));
  for jj = 1:length(fns)
    if ~isstruct(getfield(in,fns{jj}))
      if max(isnan(getfield(in,fns{jj})))
        error(['Nan detected:in.',fns{jj}])
      end
    end
  end
end

return