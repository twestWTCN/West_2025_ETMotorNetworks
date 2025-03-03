function computeFreqSpectraBank(R, sublist)
% Compute frequency spectra for the given subject list
% Input: R - structure containing paths and subject details
%        sublist - list of subjects to process

set(0, 'defaultAxesFontSize', 18);
subList = 1:numel(sublist); % Use subset if desired
prtList = 2:5; % Part list
condList = 1:4; % Condition list
parList = makeParList(subList, prtList, condList);

% Loop through virtual channels and resolution settings
for vc = 0:1
    for hr = 0:1
        parfor ip = 1:size(parList, 1)
            % Define names for different resolutions and trial types
            resname = {'lowRes', 'highRes'};
            keeptrlname = {'avtrl', 'keeptrl'};
            vcname = {'', '_VC'};
            highres = hr;

            % Set up retrieval details for data
            retdat = struct();
            retdat.sub = sublist{parList(ip, 1)};
            retdat.block = 3; % Task block
            retdat.cond = parList(ip, 3); % Condition
            retdat.part = parList(ip, 2); % Movement part
            if vc == 1
                retdat.fileapend = "[R.epoch.names{part} '_trans_pp_VC']";
            else
                retdat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']";
            end
            retdat.subfold = "['Condition' num2str(cond) '_epoched']";
            retdat.ftflag = 1;
            retdat.keeptrialinfo = 'yes';
            ftdata_cat = retrieveData(R, retdat);

            % Handle position information in virtual channels
            if isfield(ftdata_cat, 'pos')
                posflag = 1;
                postmp = ftdata_cat.pos;
                ftdata_cat = rmfield(ftdata_cat, 'pos');
            else
                posflag = 0;
                postmp = [];
            end

            % Compute spectra using multitaper method
            cfg = [];
            cfg.method = 'mtmconvol';
            cfg.taper = 'dpss';
            cfg.output = 'pow';
            if highres == 1
                cfg.toi = -4:0.05:4;
                cfg.foi = 2:0.5:98;
            else
                cfg.toi = -4:0.1:3;
                cfg.foi = 2:0.5:48;
            end
            cfg.t_ftimwin = 4 ./ cfg.foi;
            cfg.tapsmofrq = 0.4 * cfg.foi;
            cfg.keeptrials = 'yes';
            cfg.pad = 10;
            freq = ft_freqanalysis(cfg, ftdata_cat);

            % Reassign position information to frequency structure if present
            if posflag
                freq.pos = postmp;
            end

            % Add processing history if available
            if isfield(ftdata_cat.hdr, 'history')
                freq = addHistoryField(freq, 'FreqBankCompute', ftdata_cat.hdr.history);
            end

            % Save frequency data while keeping trials
            fileappend = [R.epoch.names{parList(ip, 2)} '_' resname{highres + 1} '_' keeptrlname{2} vcname{vc + 1} '_trans_freq'];
            saveExpData(R, sublist{parList(ip, 1)}, 'TaskFreq', [], fileappend, freq, ['Condition' num2str(parList(ip, 3)) '_epoched']);

            % Collapse trials
            if posflag
                freq = rmfield(freq, 'pos');
            end
            freqcol = ft_freqdescriptives([], freq);
            if posflag
                freqcol.pos = postmp;
            end

            % Save collapsed frequency data
            fileappend = [R.epoch.names{parList(ip, 2)} '_' resname{highres + 1} '_' keeptrlname{1} vcname{vc + 1} '_trans_freq'];
            saveExpData(R, sublist{parList(ip, 1)}, 'TaskFreq', [], fileappend, freqcol, ['Condition' num2str(parList(ip, 3)) '_epoched']);

            disp(['Progress: Completed ' num2str(ip) ' out of ' num2str(size(parList, 1))]);
        end
    end
end

end

% Helper function to generate parameter list
function parList = makeParList(subList, prtList, condList)
i = 0;
parList = [];
for sub = subList
    for part = prtList
        for cond = condList
            i = i + 1;
            parList(i, :) = [sub, part, cond];
        end
    end
end
end