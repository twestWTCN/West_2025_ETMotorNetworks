function [ftdata_cat,pathlist] = retrieveData(R,retdat)
if ~isfield(retdat,'keeptrialinfo')
    retdat.keeptrialinfo = 'yes';
end
if ~isfield(R,'export')
    R.export.pathonly = 0;
end


freqflag = 0;
ip = 0;
for sub = retdat.sub
    for block = retdat.block % Block of the Experiment
        switch R.import.blockpart{block}
            case 'Task'
                if strcmp(retdat.fileapend{1}(end-5:end-2),'freq')
                    fold = 'TaskFreq'; % this is frequncy domain data
                    freqflag = 1; % flag for freq data
                else
                    fold = 'Task';
                end

                for cond = retdat.cond % Condion type
                    for part = retdat.part % part of task
                        ip = ip + 1;
                        %                         fileappend = [R.epoch.names{part} '_block_pp_vc'];
                        fileappend = eval(retdat.fileapend);
                        [ftdata_rep{ip},pathlist{ip}] = loadExpData(R,sub{1},fold,[],fileappend,eval(retdat.subfold));
                    end
                end
            otherwise
                ip = ip + 1;
                %             fileappend = 'pp_vc';
                fileappend = eval(retdat.fileapend);
                [ftdata_rep{ip},pathlist{ip}] = loadExpData(R,sub{1},R.import.blockpart{block},[],fileappend,retdat.subfold);
        end
    end
end

if R.export.pathonly == 1
    ftdata_cat = [];
    return
end

if retdat.ftflag
    repN = ip;
    % Ensure same channels (or error!)
    list = [];
    for i = 1:numel(ftdata_rep)
        if i == 1
            list  =  ftdata_rep{i}.label;
        else
            list = intersect(list, ftdata_rep{i}.label);
        end
    end
    % select common channels
    cfg = [];
    cfg.channel = list;
    ip = 0;
    postmp = {}; posflag = 0;
    for i = 1:numel(ftdata_rep)
        if isfield(ftdata_rep,'trial')
            if isempty(ftdata_rep{i}.trial)
                warning('You are merging data that does not have valid trials!');
            end
        else

            if isfield(ftdata_rep{i},'time')
                xt = size(ftdata_rep{i}.time,2);
            else
                xt = 1; % hack for non time domain data
            end

            if xt>0
                ip = ip+1;
                if isfield(ftdata_rep{i},'pos')
                    % temporarily remove the pos field as messes with FT
                    postmp{ip} = ftdata_rep{i}.pos;
                    ftdata_rep{i} = rmfield(ftdata_rep{i},'pos');
                    posflag = 1;
                end
                ftdata_rep{ip} = ft_selectdata(cfg,ftdata_rep{i});
            else
                warning('Data has no valid trials!');
            end
        end
    end

    % remove trial info
    if strcmp(retdat.keeptrialinfo,'no')
        for i = 1:numel(ftdata_rep)
            ftdata_rep{i} = rmfield(ftdata_rep{i},'trialinfo');
        end
    else
        if any(cellfun(@(x) isfield(x,'trialinfo'),ftdata_rep));
            % ensure all trialinfo lengths are the same
            trialInfoN = cellfun(@(x) size(x.trialinfo,2),ftdata_rep);
            if numel(unique(trialInfoN))>1
                for i = 1:numel(ftdata_rep)
                    [nRep,nPar] = size(ftdata_rep{i}.trialinfo);
                    ftdata_rep{i}.trialinfo = [ftdata_rep{i}.trialinfo, nan(nRep,max(trialInfoN)-nPar)];
                end
            end
        end
    end

    if freqflag == 0
        ftdata_cat = appendFTData(ftdata_rep,repN);
    else
        cfg = [];
        cfg.appenddim = 'rpt';
        cfg.tolerance = 1e-3;
        ftdata_cat = ft_appendfreq(cfg,ftdata_rep{:});
    end

    % For Pos info
    if posflag == 1
        if numel(unique(cellfun(@length,postmp))) > 1
            warning('Your Position fields have different lengths!!!')
        end
        ftdata_cat.pos = postmp{1}; % fill in with first instance
    end
else
    ftdata_cat = [ftdata_rep{:}];
end

