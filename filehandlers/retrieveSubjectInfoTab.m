function T = retrieveSubjectInfoTab(R,sublist)
pathn = [R.path.datapath '\' R.path.expname '\Group\SubjectDetails.xlsx'];
T = readtable(pathn,'NumHeaderLines',0,'VariableNamingRule','preserve');
[~,IA,IB] = intersect(sublist',T.Subject,'stable');
T = T(IB,:);
