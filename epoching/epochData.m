function [] = epochData(R)
% epochData: Epochs the steady state blocks (rest/posture) into regular-sized chunks.
% The task block is further subdivided into three blocks: Motor preparation, motor execution, and hold.
% Parameters:
%   R - Configuration structure containing various settings and parameters for the experiment.
% This function also handles subdivisions of blocks as well as transitions between task states.

% Define which method to use for epoching
R.epoch.manfresh = 3; % 3: Using pre-trained CNN for automated epoching

% Load pre-trained CNN networks if method 3 is chosen
if R.epoch.manfresh == 3
    netHold = loadExpData(R, 'Group', 'AutoMoveLocking', R.epoch.names{5}, 'MarkerNet', []);
    netExec = loadExpData(R, 'Group', 'AutoMoveLocking', R.epoch.names{4}, 'MarkerNet', []);
    netQual = loadExpData(R, 'Group', 'AutoMoveLocking', 'QualityNet', 'QualityNet', []);
end

blockpart = R.epoch.blockpart;
subI =0;
for sub = R.epoch.subsel
    subI = subI + 1;
    for block =1:numel(blockpart) % Block of the Experiment
        switch blockpart{block}
            case {'Noise','Rest','Posture'}
                [ftdata,~,prsphist] = loadExpData(R,sub{1},blockpart{block},[],'pp','pp');

                %bugfix for EEG to eeg
                ftdata.hdr.chantype(strcmp(ftdata.hdr.chantype,'accelerometer')) = {'Accelerometer'};

                ftdata = cleanSteadyStateByAcc(ftdata,95); % this function makes sure only periods of stable ACC are included
                cfg = [];
                cfg.length = R.epoch.length;
                cfg.keepsampleinfo = 'no';
                ftdata = ft_redefinetrial(cfg,ftdata); % Split up
                ftdata = removeBadTrials(ftdata);
                ftdata = repairBadChannel(R,ftdata);

                ftdata = addHistoryField(ftdata,'Epoched',prsphist); % add history field
                saveExpData(R,sub{1},blockpart{block},[],'pp',ftdata,'epoched');
            case 'Task'
                % Load Each Condition
                [ftdata,~,prsphist] = loadExpData(R,sub{1},blockpart{block},[],'pp','pp');

                % Bugfix for EEG header chantype from 'accelerometer' to 'Accelerometer'
                ftdata.hdr.chantype(strcmp(ftdata.hdr.chantype,'accelerometer')) = {'Accelerometer'};

                % Find Coder Channel Index
                coderInd = find(strcmp(ftdata.label,'Coder')); % Use the local coder as it is more accurate

                %% Concatenate All Repetitions into One Large ftdata Structure
                ftdata_tmp = ftdata;
                ftdata_tmp.trial = {[ftdata.trial{:}]};
                ftdata_tmp.time = {linspace(0, size(ftdata_tmp.trial{1}, 2) / ftdata_tmp.fsample, size(ftdata_tmp.trial{1}, 2))};
                if isfield(ftdata_tmp, 'trialinfo')
                    ftdata_tmp = rmfield(ftdata_tmp, 'trialinfo');
                end

                %% Replace NaNs with Zeros (necessary to avoid deletion by FieldTrip)
                X = ftdata_tmp.trial{1};
                X(isnan(X)) = 0;
                ftdata_tmp.trial{1} = X;

                %% Extract Event Definitions from the 'Coder' Channel
                coder = ftdata_tmp.trial{1}(coderInd, :);

                % Define Task Codes
                % KEY CODES
                %   1: pretrial rest               = 0.2v
                %   2: pretrial posture            = 0.5v
                %   3: movement cue (black arrows) = 1.5v
                %   4: go signal (green arrows)    = 2.0v
                %   5: hold (collision)            = 2.5v
                %   6: successful hold (pop)       = 3v

                switch R.import.site
                    case 'UCL'
                        levels = [0.2, 0.5, 1.5, 2.0, 2.5]'; % rest, post, prep, exec, hold
                    case 'DL'
                        levels = [0.002, 0.005, 0.015, 0.020, 0.025]'; % rest, post, prep, exec, hold
                end

                for cond = R.import.condsel
                    % Subdivide Task
                    trialdetails = {{['Condition ' num2str(cond)]}, {sub{1}}};
                    minblocksize = [0.5, 0.5, 0.1, 0, 0.5] * ftdata_tmp.fsample;
                    maxblocksize = [3.1, 3.1, inf, inf, inf] * ftdata_tmp.fsample;
                    maxbreaksize = [1, 1, 1, 7, 7] * ftdata_tmp.fsample;
                    [decoder, SMRtrialdef] = codetimings(levels, coder, trialdetails, 'taskSplit', minblocksize, maxbreaksize, [], maxblocksize);

                    % Load Unity Data (trial information)
                    UNData = loadExpData(R, sub{1}, 'Unity', [], '', 'table');
                    UNConds = table2cell(UNData(:, [1, 3, 5]));
                    epochSanityCheck(R, SMRtrialdef, UNData, ftdata); % Sanity check for epoching
                    trialdetails = {{['Condition ' num2str(cond)]}, {sub{1}}};

                    % Condition Split
                    switch R.import.site
                        case 'UCL'
                            levelsC = [3.5, 4, 4.5, 5]'; % Condition 1, 2, 3, 4
                        case 'DL'
                            levelsC = [0.035, 0.040, 0.045, 0.050]'; % Condition 1, 2, 3, 4
                    end
                    [decoder, condtrialdef] = codetimings(levelsC, coder, trialdetails, 'condSplit', [0.5, 0.5, 0.5, 0.5] * ftdata_tmp.fsample, [1, 1, 1, 1] * ftdata_tmp.fsample, UNConds);

                    % Merge Condition Information with Task Data
                    SMRtrialdef = selectTrialsByCondition(SMRtrialdef, condtrialdef{cond, 2}, condtrialdef{cond, 5});
                    for part = 1:5
                        % Construct Trial Data
                        if part < 6
                            trialInfo{part} = [cell2mat(SMRtrialdef{part, 5}(1, :)); cell2mat(SMRtrialdef{part, 5}(2, :)); grp2idx([SMRtrialdef{part, 5}{3, :}])' - 2]';
                        else
                            % For the full sequence data - everything is defined from the start of the rest trials.
                            % Hold is not always met
                            trialInfo{part} = [cell2mat(SMRtrialdef{1, 5}(1, :)); cell2mat(SMRtrialdef{1, 5}(2, :)); grp2idx([SMRtrialdef{1, 5}{3, :}])' - 2]';
                        end

                        baseMed = [];
                        for prt = 1:5
                            baseMed(:, prt) = getTrlAverageBaseline(ftdata_tmp, SMRtrialdef{prt, 2}, R.epoch.accBadFlag{1}(subI));
                            scatter(prt, baseMed(:, prt), 'MarkerFaceColor', R.plot.cmap.lplots(prt, :), 'MarkerEdgeColor', 'none'); hold on;
                        end
                        clear prt

                        switch R.import.site
                            case 'UCL'
                                vecthresh = 0.25;
                            case 'DL'
                                vecthresh = 0.25;
                        end

                        %% Get Transition Data
                        epcFunPost = @convertTimeLockToPost_AUTO;

                        cfg = [];
                        tranSamp = fix((R.epoch.trans{part} / 1000) * ftdata_tmp.fsample);
                        blockSamp = fix((R.epoch.block{part} / 1000) * ftdata_tmp.fsample);
                        if part == 1 % Rest
                            searchwin = R.epoch.searchwin{1};
                            searchwin = fix((searchwin / 1000) * ftdata_tmp.fsample);
                            cfg.trlBlock = [SMRtrialdef{1, 2}(1, :) + searchwin(1); SMRtrialdef{1, 2}(2, :) + searchwin(2); repmat(searchwin(1), size(SMRtrialdef{1, 2}(1, :))); trialInfo{1}']';
                            cfg.trlTrans = cfg.trlBlock;
                        elseif part == 2 % Rest -> Prep (Hands Up)
                            searchwin = R.epoch.searchwin{2};
                            searchwin = fix((searchwin / 1000) * ftdata_tmp.fsample);
                            cfg = [];
                            cfg.trlTrans = [SMRtrialdef{1, 2}(2, :) + searchwin(1); SMRtrialdef{1, 2}(2, :) + searchwin(2); repmat(searchwin(1), size(SMRtrialdef{1, 2}(1, :))); trialInfo{1}']';
                            reactRange = [-0.5, 4] * 1000;
                            close all;
                            [cfg.trlTrans, reacttime] = epcFunPost(ftdata_tmp, cfg.trlTrans, tranSamp, 1, R.epoch.accBadFlag{1}(subI), searchwin, vecthresh, reactRange, baseMed(:, 1));
                            cfg.trlBlock = [SMRtrialdef{1, 2}(2, :) + searchwin(1); SMRtrialdef{1, 2}(2, :) + searchwin(2); repmat(searchwin(1), size(SMRtrialdef{1, 2}(1, :))); trialInfo{1}']';
                            reactbank{part, cond} = reacttime;
                            clf;
                        elseif part == 3 % Posture -> Prep (Arrows Appear)
                            cfg = [];
                            cfg.trlTrans = [SMRtrialdef{2, 2}(2, :) + tranSamp(1); SMRtrialdef{2, 2}(2, :) + tranSamp(2); repmat(tranSamp(1), size(SMRtrialdef{2, 2}(1, :))); trialInfo{2}']';
                            cfg.trlBlock = [SMRtrialdef{2, 2}(2, :) + blockSamp(1); SMRtrialdef{2, 2}(2, :) + blockSamp(2); repmat(blockSamp(1), size(SMRtrialdef{2, 2}(1, :))); trialInfo{2}']';
                            reactbank{part, cond} = ones(1, size(cfg.trlTrans, 1));
                            clf;
                        elseif part == 4 % Prep -> Exec
                            searchwin = R.epoch.searchwin{4};
                            searchwin = fix((searchwin / 1000) * ftdata_tmp.fsample);
                            if R.epoch.manfresh == 1
                                cfg = [];
                                cfg.trlTrans = [SMRtrialdef{4, 2}(1, :) + searchwin(1); SMRtrialdef{4, 2}(1, :) + searchwin(2); repmat(searchwin(1), size(SMRtrialdef{4, 2}(1, :))); trialInfo{4}']';
                                [exec_trl, exec_hold, reacttime_exec_prephold, testdata] = convertTimeLockToExec_MAN(ftdata_tmp, cfg.trlTrans, tranSamp, 1, R.epoch.accBadFlag{1}(subI), searchwin);
                                cfg.trlTrans = exec_trl;
                            elseif R.epoch.manfresh == 2
                                try
                                    fileappend = [R.epoch.names{part} '_trlTransDef_V3'];
                                    cfg.trlTrans = loadExpData(R, sub{1}, 'Task', [], fileappend, ['Condition' num2str(cond) '_epoched']);
                                    reacttime_exec_prephold = {ones(1, size(cfg.trlTrans, 1)), ones(1, size(cfg.trlTrans, 1))};
                                catch
                                    warning('Cant load prelabelled trial conditions!');
                                    R.epoch.manfresh = 1;
                                end
                            elseif R.epoch.manfresh == 3
                                cfg = [];
                                cfg.trlTrans = [SMRtrialdef{4, 2}(1, :) + searchwin(1); SMRtrialdef{4, 2}(1, :) + searchwin(2); repmat(searchwin(1), size(SMRtrialdef{4, 2}(1, :))); trialInfo{4}']';
                                [exec_trl, exec_hold, reacttime_exec_prephold] = convertTimeLockToMarker_CNN(ftdata_tmp, cfg.trlTrans, tranSamp, 1, R.epoch.accBadFlag{1}(subI), searchwin, netExec, netHold, netQual);
                                cfg.trlTrans = exec_trl;
                            end
                            cfg.trlBlock = [SMRtrialdef{4, 2}(1, :) + searchwin(1); SMRtrialdef{4, 2}(1, :) + searchwin(2); repmat(searchwin(1), size(SMRtrialdef{4, 2}(1, :))); trialInfo{4}']';
                            reactbank{part, cond} = reacttime_exec_prephold{1};
                            clf;
                        elseif part == 5 % Exec -> Hold
                            searchwin = R.epoch.searchwin{5};
                            if R.epoch.manfresh == 1
                                cfg.trlTrans = exec_hold;
                            elseif R.epoch.manfresh == 2
                                fileappend = [R.epoch.names{part} '_trlTransDef_V3'];
                                cfg.trlTrans = loadExpData(R, sub{1}, 'Task', [], fileappend, ['Condition' num2str(cond) '_epoched']);
                            elseif R.epoch.manfresh == 3
                                cfg.trlTrans = exec_hold;
                            end
                            reactbank{part, cond} = reacttime_exec_prephold{2};
                            cfg.trlBlock = [SMRtrialdef{5, 2}(1, :) + tranSamp(1); SMRtrialdef{5, 2}(1, :) + tranSamp(2); repmat(tranSamp(1), size(SMRtrialdef{5, 2}(1, :))); trialInfo{5}']';
                            clf;
                        end
                        trlBlockSave = cfg.trlBlock;
                        if part > 1 && part < 6
                            trlTransSave = cfg.trlTrans;
                            cfg.trl = cfg.trlTrans;
                            ftdata_trans = ft_redefinetrial(cfg, ftdata_tmp);
                            ftdata_trans = repairBadChannel(R, ftdata_trans);
                            if part < 6
                                checkEpochAcc(ftdata_trans, R.epoch.accBadFlag{1}(subI));
                            else
                                checkEpochAccFullSeq(ftdata_trans, R.epoch.accBadFlag{1}(subI));
                            end
                        end
                        cfg.trl = cfg.trlBlock;
                        ftdata_block = ft_redefinetrial(cfg, ftdata_tmp);
                        [ftdata_block] = removeBadTrials(ftdata_block);
                        ftdata_block = repairBadChannel(R, ftdata_block);

                        % Save Block Data
                        ftdata_block = addHistoryField(ftdata_block, 'Epoched', prsphist);
                        ft_nancheck(ftdata_block);
                        fileappend = [R.epoch.names{part} '_blockMot_ppV3'];
                        saveExpData(R, sub{1}, 'Task', [], fileappend, ftdata_block, ['Condition' num2str(cond) '_epoched']);
                        fileappend = [R.epoch.names{part} '_trlBlockDef_V3'];
                        saveExpData(R, sub{1}, 'Task', [], fileappend, trlBlockSave, ['Condition' num2str(cond) '_epoched']);

                        % Save Transition data
                        if part > 1 && part < 6
                            ftdata_trans = addHistoryField(ftdata_trans, 'Epoched', prsphist);
                            ft_nancheck(ftdata_trans);
                            fileappend = [R.epoch.names{part} '_transMot_ppV3'];
                            checkConsistentDataLength(ftdata_trans, 0.1);
                            saveExpData(R, sub{1}, 'Task', [], fileappend, ftdata_trans, ['Condition' num2str(cond) '_epoched']);
                            fileappend = [R.epoch.names{part} '_trlTransDef_V3'];
                            saveExpData(R, sub{1}, 'Task', [], fileappend, trlTransSave, ['Condition' num2str(cond) '_epoched']);
                        end
                    end
                end
        end % Block Loop
    end % Subject Loop
end
end
function [ftdata, reactbank] = removeBadTrials(ftdata, reactbank)
% REMOVEBADTRIALS removes trials with bad (NaN) values from the dataset.
% This function identifies trials with NaN values and removes them from both the ftdata and reactbank.

badtr = [];
for tr = 1:numel(ftdata.trial)
    X = ftdata.trial{tr};
    badtr(tr) = any(isnan(X(:))); % Identify trials containing NaN values
end

fprintf('Removing %.f trials containing visually marked artefacts \n', sum(badtr));

% Select good trials only
cfg = [];
cfg.trials = find(~badtr);
ftdata = ft_selectdata(cfg, ftdata);

% Update reactbank if provided
if nargin > 1
    reactbank = reactbank(~badtr);
end

end

function epochSanityCheck(R, SMRtrialdef, UNData, ftdata)
% EPOCHSANITYCHECK verifies if epoch definitions match the unity game data.
% This function checks that the number of epochs and their lengths are consistent with expected values.

% Check if number of epochs matches the unity game data
if size(SMRtrialdef{1, 2}, 2) ~= size(UNData, 1)
    warning('Your epoch definitions do not match the unity game data! Will try and fix in cond epoching...');
end

% Check lengths of rest and posture trials
sectionLength = [3 3]; % Length of rest and posture trials (in seconds)
for seg = 1:2 % Loop over rest and posture trials
    if any((diff(SMRtrialdef{seg, 2}) - sectionLength(seg) * ftdata.fsample) > (0.1 * ftdata.fsample))
        warning(['Epoch lengths for ' R.epoch.names{seg} ' are of unexpected length!']);
    end
end

end

function ftdata = repairBadChannel(R, ftdata)
% REPAIRBADCHANNEL identifies and repairs bad channels in EEG/MEG data.
% This function searches for channels containing NaN values, identifies them as bad, and attempts to repair them using data from neighboring channels.

% Find channels with NaN values (indicating bad channels)
badch = find(sum(isnan(ftdata.trial{1}), 2));

% If bad channels are found, proceed with repair
if ~isempty(badch)
    cfg = [];
    cfg.badchannel = ftdata.label(badch); % List of bad channels to repair
    cfg.neighbours = R.neighbours; % Neighbors definition for interpolation
    ftdata = ft_channelrepair(cfg, ftdata); % Repair bad channels using neighboring data
else
    % No bad channels found, return data as is
    ftdata = ftdata;
end

end
