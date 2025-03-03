function [] = artefactRejectData(R)
% ARTEFACTREJECTDATA processes and removes artefacts from experimental data.
% This function loads data, applies artefact rejection, and saves cleaned data.
% It operates over various conditions, subjects, and blocks of experimental data.

subI = 0;
for sub = R.artRej.subsel
    subI = subI + 1;
    for block = 1:numel(R.artRej.blockpart) % block loop (rest/posture/task)
        switch R.artRej.blockpart{block}
            case 'Task'
                prun = 5;
                condlist = R.import.condsel;
            otherwise
                prun = 1;
                condlist = 1;
        end

        %% Main Loop to Collect Data
        catI = 0; ftdata_tmp = {}; catDAT = {}; filelist = [];
        for cond = condlist
            for part = 1:prun % task part loop
                if part > 1 && part < 6
                    ftrun = 2;
                else
                    ftrun = 1; % Motor prep doesn't have a transition
                end
                for ft = 1:ftrun % filetypes block, timelocked
                    catI = catI + 1;
                    R.path.tracker = [sub cond block part ft];
                    filelist(catI, :) = [cond part ft];
                    switch R.artRej.blockpart{block}
                        case 'Task'
                            if ft == 1
                                fileappend = [R.epoch.names{part} '_blockMot_ppV3'];
                                ftdata_tmp{catI} = loadExpData(R, sub{1}, 'Task', [], fileappend, ['Condition' num2str(cond) '_epoched']);
                            elseif ft == 2
                                fileappend = [R.epoch.names{part} '_transMot_ppV3'];
                                ftdata_tmp{catI} = loadExpData(R, sub{1}, 'Task', [], fileappend, ['Condition' num2str(cond) '_epoched']);
                            end
                            catDAT{catI} = repmat(filelist(catI, :), numel(ftdata_tmp{catI}.trial), 1);
                            if isfield(ftdata_tmp{catI}, 'trialinfo')
                                trialinfoCat{catI} = ftdata_tmp{catI}.trialinfo;
                            else
                                disp('No trialinfo available!')
                                trialinfoCat{catI} = repmat([cond 0 0], numel(ftdata_tmp{catI}.trial), 1);
                            end
                        otherwise
                            fileappend = 'pp';
                            ftdata_tmp{catI} = loadExpData(R, sub{1}, R.artRej.blockpart{block}, [], fileappend, 'epoched');
                            catDAT{catI} = repmat(filelist(catI, :), numel(ftdata_tmp{catI}.trial), 1);
                            trialinfoCat{catI} = repmat([cond 0 0], numel(ftdata_tmp{catI}.trial), 1);
                    end
                end
            end
        end
        
        % Temporarily Remove Trial Info (replace it later)
        if strcmp(R.artRej.blockpart{block}, 'Task')
            for i = 1:numel(ftdata_tmp)
                if isfield(ftdata_tmp{i}, 'trialinfo')
                    ftdata_tmp{i} = rmfield(ftdata_tmp{i}, 'trialinfo');
                end
            end
        end

        % Append data together
        ftdata_cat = ft_appenddata([], ftdata_tmp{:});
        ftdata_cat.hdr = ftdata_tmp{1}.hdr;
        
        % Check Headers
        checkHeader(ftdata_cat);
        catdat = vertcat(catDAT{:});

        % Ensure all trialinfos are the same dimension to concatenate
        maxTr = max(cellfun(@(x) size(x, 2), trialinfoCat));
        trialinfoCat = cellfun(@(x) fillTrialInfo(x, maxTr), trialinfoCat, 'UniformOutput', 0);
        trialinfodat = vertcat(trialinfoCat{:});

        %% Run the Artefact Rejection
        switch R.import.site
            case 'UCL'
                [ftdata_cat, badTrials] = artefact_rejection(R, ftdata_cat, {'tangental', 'radial', 'axial'});
            case 'DL'
                [ftdata_cat, badTrials] = artefact_rejection(R, ftdata_cat, {'EEG'});
        end
        
        % Remove bad trials from data
        catdat(badTrials, :) = [];
        trialinfodat(badTrials, :) = [];

        %% Repartition Data
        for cond = condlist
            for part = 1:prun % task part loop
                if part > 1 && part < 6
                    ftrun = 2;
                else
                    ftrun = 1; % Motor prep doesn't have a transition
                end
                for ft = 1:ftrun % filetypes block, timelocked
                    trSel = find(all(catdat == [cond part ft], 2));
                    if numel(trSel) < 16
                        warning('Less than 16 trials remain after AR!')
                    end

                    cfg = [];
                    cfg.trials = trSel;
                    ftdata = ft_selectdata(cfg, ftdata_cat);
                    ftdata.trialinfo = trialinfodat(trSel, :);
                    
                    % Check Headers
                    checkHeader(ftdata);
                    
                    if part ~= 6 % Full reach is different length
                        ftdata = checkConsistentDataLength(ftdata, 0.2);
                    end
                    
                    %% Save Data
                    ftdata = addHistoryField(ftdata, 'ArtRej');
                    switch R.artRej.blockpart{block}
                        case 'Task'
                            if ft == 1
                                fileappend = [R.epoch.names{part} '_block_pp_arV3'];
                                saveExpData(R, sub{1}, 'Task', [], fileappend, ftdata, ['Condition' num2str(cond) '_epoched']);
                            elseif ft == 2
                                fileappend = [R.epoch.names{part} '_trans_pp_arV3'];
                                saveExpData(R, sub{1}, 'Task', [], fileappend, ftdata, ['Condition' num2str(cond) '_epoched']);
                            end
                        otherwise
                            fileappend = 'pp_arV3';
                            saveExpData(R, sub{1}, R.artRej.blockpart{block}, [], fileappend, ftdata, 'epoched');
                    end
                    disp([sub cond block part]);
                end % epoch type loop
            end % task part loop
        end % cond loop
    end % block loop
end % Sub loop
end
