function [object,pathload,prsphist] = loadExpData(R,sub,blockpart,condname,ext,subfold)
if ~isfield(R,'export')
    R.export.pathonly = 0;
end
% This function will load in data for the analysis
% It has the path format - root\experimentName\session\
fileset = [sub '_' blockpart];
if ~isempty(condname)
    fileset = [fileset '_' condname];
end

root = [R.path.datapath '\' R.path.expname '\' sub ...
    '\Fieldtrip\' fileset '\' subfold '\'];
% mkdir(root)
pathload = [root fileset '_' ext];
if R.export.pathonly == 1
    pathload = [pathload '.mat'];
    object = [];
    return
end
object = load(pathload,'target');
object = object.target;

if isfield(object,'hdr')
    if isfield(object.hdr,'history')
        prsphist = object.hdr.history; % keep history persistant
    end
end