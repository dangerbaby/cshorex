[fList,pList] = matlab.codetools.requiredFilesAndProducts('explore_prediction.m');
for i = 1:length(fList)
  if ~isempty(strfind(fList{i},'matlab'))
    copyfile(fList{i},'./mfiles/')
  end
end