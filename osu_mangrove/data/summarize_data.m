%https://www.designsafe-ci.org/data/browser/public/designsafe.storage.published/PRJ-2143/NHERIMangroves/data/inter
clear all
basename = 'HighDensity_h270_hv182_NoWall'
%basename = 'Baseline_h270_hv185_NoWall'
dnames = dir(['./',basename,'/T*']);

%dname = './HighDensity_h270_hv182_NoWall/Trial08/';
%dname = './HighDensity_h270_hv182_NoWall/Trial09/';

for j = 1:length(dnames)

  dname = ['./',basename,'/',dnames(j).name,'/']
  %pressure 
  fnames = dir([dname,'press*']);
  %fnames.name
  for i = 1:length(fnames)
    %for i = 1
    id = str2num(fnames(i).name(6:find(fnames(i).name=='.')-1));
    dat.press(id).name = fnames(i).name;
    fid = fopen([dname,fnames(i).name]);
    dum = 'abc';
    while ~contains(dum,'[Data]')
      dum = fgetl(fid);
      if contains(dum,'StillWaterDepth')
        dat.press(id).swd = str2num(dum(19:end)); 
      end
      if contains(dum,'X:')
        dat.press(id).x = str2num(dum(5:end)); 
      end
      if contains(dum,'Z:')
        dat.press(id).z = str2num(dum(5:end)); 
      end
      if contains(dum,'[Data]')
        dat.press(id).press = cell2mat(textscan(fid,'%f','Delimiter', '\n'));
      end
    end
    fclose(fid);
  end

  %wg
  fnames = dir([dname,'wg*']);
  for i = 1:length(fnames)
    %for i = 1
    id = str2num(fnames(i).name(3:find(fnames(i).name=='.')-1));
    dat.wg(id).name = fnames(i).name;
    fid = fopen([dname,fnames(i).name]);
    dum = 'abc';
    while ~contains(dum,'[Data]')
      dum = fgetl(fid);
      if contains(dum,'X:')
        dat.wg(id).x = str2num(dum(5:end)); 
      end
      if contains(dum,'[Data]')
        dat.wg(id).eta = cell2mat(textscan(fid,'%f','Delimiter', '\n'));
      end
    end
    fclose(fid);
  end

  %uswg
  fnames = dir([dname,'uswg*']);
  for i = 1:length(fnames)
    %for i = 1
    id = str2num(fnames(i).name(5:find(fnames(i).name=='.')-1));
    dat.uswg(id).name = fnames(i).name;
    fid = fopen([dname,fnames(i).name]);
    dum = 'abc';
    while ~contains(dum,'[Data]')
      dum = fgetl(fid);
      if contains(dum,'X:')
        dat.uswg(id).x = str2num(dum(5:end)); 
      end
      if contains(dum,'[Data]')
        dat.uswg(id).eta = cell2mat(textscan(fid,'%f','Delimiter', '\n'));
      end
    end
    fclose(fid);
  end

  %u
  clear fnames;
  fnames2 = dir([dname,'u*.txt']);
  cnt = 0;
  for i = 1:length(fnames2)
    if ~contains(fnames2(i).name,'wg')
      cnt = cnt+1;
      fnames(cnt)=fnames2(i);
    end
  end
  for i = 1:length(fnames)
    %for i = 1
    id = str2num(fnames(i).name(2:find(fnames(i).name=='.')-1));
    dat.u(id).name = fnames(i).name;
    fid = fopen([dname,fnames(i).name]);
    dum = 'abc';
    while ~contains(dum,'[Data]')
      dum = fgetl(fid);
      if contains(dum,'X:')
        dat.u(id).x = str2num(dum(5:end)); 
      end
      if contains(dum,'Z:')
        dat.u(id).z = str2num(dum(5:end)); 
      end
      if contains(dum,'[Data]')
        dat.u(id).u = cell2mat(textscan(fid,'%f','Delimiter', '\n'));
      end
    end
    fclose(fid);
  end

  %w
  clear fnames;
  fnames2 = dir([dname,'w*.txt']);
  cnt = 0;
  for i = 1:length(fnames2)
    if ~contains(fnames2(i).name,'wg')&~contains(fnames2(i).name,'wm')
      cnt = cnt+1;
      fnames(cnt)=fnames2(i);
    end
  end
  for i = 1:length(fnames)
    %for i = 1
    id = str2num(fnames(i).name(2:find(fnames(i).name=='.')-1));
    dat.w(id).name = fnames(i).name;
    fid = fopen([dname,fnames(i).name]);
    dum = 'abc';
    while ~contains(dum,'[Data]')
      dum = fgetl(fid);
      if contains(dum,'X:')
        dat.w(id).x = str2num(dum(5:end)); 
      end
      if contains(dum,'Z:')
        dat.w(id).z = str2num(dum(5:end)); 
      end
      if contains(dum,'[Data]')
        dat.w(id).w = cell2mat(textscan(fid,'%f','Delimiter', '\n'));
      end
    end
    fclose(fid);
  end

  save([dname,'summary.mat'],'dat')
end
