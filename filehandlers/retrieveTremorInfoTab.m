function T = retrieveTremorInfoTab(R,sublist)
if strcmp(sublist{1}(1:2),'OP')
    pathn = [R.path.datapath '\' R.path.expname '\Group\TremorAssessmentData.xlsx'];
else
    pathn = [R.path.datapath '\' R.path.expname '\Group\TremorAssessmentData.xlsx'];
end
    T = readtable(pathn,'NumHeaderLines',0,'VariableNamingRule','preserve');
    [~,IA,IB] = intersect(sublist',T.Name,'stable');
    T = T(IB,:);
