function pathsave = saveExpData(R,sub,blockpart,condname,ext,target,subfold)
if isfield(target,'cfg')
    target = rmfield(target,'cfg');
end
if isfield(target,'hdr')
    if isfield(target.hdr,'label')
        checkHeader(target)
    end
end
fileset = [sub '_' blockpart];
if ~isempty(condname)
    fileset = [fileset '_' condname];
end

root = [R.path.datapath '\' R.path.expname '\' sub ...
    '\Fieldtrip\' fileset '\' subfold '\'];
mkdir(root)
pathsave = [root fileset '_' ext];
if ~isempty(target); %allows to just recover path
    s = whos('target');
    if (s.bytes/1e9) < 2
        save(pathsave,'target')
    else
        warning('Using v7.3 save, this can be very slow!')
        save(pathsave,'target','-v7.3')
    end
end