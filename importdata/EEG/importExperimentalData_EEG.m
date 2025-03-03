function importExperimentalData_EEG(R)
    % IMPORTEXPERIMENTALDATA_EEG: Imports EEG, AV, and motion tracking data.
    % Converts data into FieldTrip format, aligns data sets, and saves the processed information.
    % Parameters:
    %   R - Structure containing paths and import options

    close all;
    %% Step 1: Initialize Variables
    blockpart = R.import.blockpart;
    subI = 0;
    for sub = R.import.subsel
        subI = subI + 1;
        UIflag = 1;  % Flag for UI coder selection
        SPLITflag = 1;  % Flag for forcing manual data split
        ftdata_rep = struct([]);  % To concatenate repetitions

        %% Step 2: Iterate Through Repetitions
        for rep = R.import.repsel{subI}
            % Load AV Data
            AVFN = [R.path.datapath '\' R.path.expname '\' sub{1} '\StimuliPCLocal\' sub{1} '_Task_' num2str(rep) '\Data.txt'];

            if ~exist([fileparts(AVFN) '\AVData.mat'], 'file')
                % Load and convert AV Data
                importInitialStimPCRaw(AVFN);
                try
                    AVData = importLocalDataGerman(AVFN, R.import.rs_fs);
                catch
                    disp('Trying UK format');
                    AVData = importEngAVData(AVFN, R.import.rs_fs);
                end
            else
                load([fileparts(AVFN) '\AVData.mat'], 'AVData');
            end

            % Load EEG Data
            EEGFN = [R.path.datapath '\' R.path.expname '\' sub{1} '\EEG\' sub{1} 'OptiTrack_Rep' num2str(rep) '.mat'];

            if ~exist(EEGFN, 'file')
                if rep == R.import.repsel{subI}(1) && SPLITflag
                    BAFN = [R.path.datapath '\' R.path.expname '\' sub{1} '\EEG\' sub{1} 'OptiTrack.eeg'];
                    importBrainAmpEEGData(R, BAFN, 1);
                end
                EEGData = importBrainAmpEEGData(R, EEGFN, 0);
            else
                load(EEGFN, 'trialData');
                EEGData = trialData;
            end

            % Load MOTIVE Data
            TAKFN = [R.path.datapath '\' R.path.expname '\' sub{1} '\StimuliPCLocal\Motive\' sub{1} '_Rep_' num2str(rep) '.tak'];
            keycodepath = [R.path.datapath '\' R.path.expname '\' sub{1} '\StimuliPCLocal\' sub{1} '_Calibration\'];

            if exist(TAKFN, 'file')
                visSplit = R.import.TAK.visplit(subI);
                if ~exist([fileparts(TAKFN) '\TAKData_' num2str(rep) '.mat'], 'file')
                    TAKData = importMotiveData(TAKFN, R.import.rs_fs, R.path.takconvpath, keycodepath, visSplit);
                else
                    load([fileparts(TAKFN) '\TAKData_' num2str(rep) '.mat'], 'TAKData');
                end
            end

            %% Step 3: Pad Data for Uniform Length
            n_AV = numel(AVData.time{1});
            n_EEG = numel(EEGData.time{1});

            if exist('TAKData', 'var')
                n_TAK = numel(TAKData.time{1});
                datSize = [n_AV, n_EEG, n_TAK];
                timers = {AVData.time, EEGData.time, TAKData.time};
            else
                datSize = [n_AV, n_EEG];
                timers = {AVData.time, EEGData.time};
            end

            [~, dL] = max(datSize);

            % Pad AV Data
            if n_AV < datSize(dL)
                dSize = datSize(dL) - n_AV;
                AVData = padData(AVData, dSize, timers{dL}, 0);
            end

            % Pad EEG Data
            if n_EEG < datSize(dL)
                dSize = datSize(dL) - n_EEG;
                EEGData = padData(EEGData, dSize, timers{dL}, 0);
            end

            % Pad TAK Data (if applicable)
            if exist('TAKData', 'var') && n_TAK < datSize(dL)
                dSize = datSize(dL) - n_TAK;
                TAKData = padData(TAKData, dSize, timers{dL}, 0);
            end

            %% Step 4: Realign Data (Aligns AV and EEG)
            time_EEG = EEGData.time{1};
            trig_EEG = -EEGData.trial{1}(strncmp(EEGData.label, 'LabJack', 3), :);
            time_AV = AVData.time{1};
            trig_AV = AVData.trial{1}(strncmp(AVData.label, 'Coder', 3), :);
            trig_AV(trig_AV > 0.1) = 0.03;

            tarInd = [find(strncmp(EEGData.label, 'LabJack', 3)), find(strncmp(AVData.label, 'Coder', 3))];
            [EEGData.trial{1}, AVData.trial{1}] = realignDataV2(time_EEG, trig_EEG, time_AV, trig_AV, R.import.rs_fs, EEGData.trial{1}, AVData.trial{1}, 1, tarInd);

            % Set Data Lengths Equal
            [EEGData, AVData] = setFTLengthToMin(EEGData, AVData);

            %% Step 5: Replace Padding with NaNs
            EEGData.trial{1}(EEGData.trial{1} == 0) = nan;
            AVData.trial{1}(AVData.trial{1} == 0) = nan;

            %% Step 6: Concatenate Data Sets
            chantypeConcat = [EEGData.chantype, AVData.chantype];
            EEGData = rmfield(EEGData, 'chantype');
            AVData = rmfield(AVData, 'chantype');

            cfg = [];
            ftdata = ft_appenddata(cfg, EEGData, AVData);
            ftdata.fsample = R.import.rs_fs;
            ftdata.hdr.chantype = chantypeConcat;
            ftdata.hdr.label = ftdata.label;

            % Plot Alignment Results
            figure;
            subplot(2, 1, 1);
            cfg = [];
            cfg.channel = {'Coder', 'LabJack'};
            seldat = ft_selectdata(cfg, ftdata);
            plot(seldat.time{1}, normalize(seldat.trial{1}, 2));
            title('AV and LabJack Synchronization');

            subplot(2, 1, 2);
            cfg = [];
            cfg.channel = {'AccX', 'xFinger1Pos'};
            seldat = ft_selectdata(cfg, ftdata);
            plot(seldat.time{1}, [1; 1] .* normalize(seldat.trial{1}, 2));
            title('Motion Tracking and Finger Position');

            %% Step 7: Subdivide Data by Condition
            coder = ftdata.trial{1}(strcmp(ftdata.label, 'Coder'), :);
            if UIflag == 1
                disp(['Base', 'Rest', 'Posture', 'Task']);
                Y = [0, 0.002, 0.005, 0.025]';
                disp(Y);
                UIflag = 0;
            end

            minblocksize = [5, 20, 10, 0.1] * ftdata.fsample;
            minbreak = [5, 5, 5, 5];
            trialdetails = [{['Rep ' num2str(rep)]}, {['Condition All']}, {sub{1}}];
            [decoder, SMRtrialdef] = codetimingsEEG(Y, coder, trialdetails, 'blockSplit', minblocksize, minbreak);

            % Create Temporary Copy for Redefining Trials
            ftdata_tmp = ftdata;
            for block = 1:3  % Rest, Posture, Task
                cfg = [];
                cfg.trl = reshape(SMRtrialdef{block, 2}(:), [], 2);
                cfg.trl = [cfg.trl, zeros(size(cfg.trl, 1), 1)];
                ftdata_rep{rep, block} = ft_redefinetrial(cfg, ftdata_tmp);
            end

            % Clear Variables for Next Iteration
            clear ftdata_tmp AVData EEGData TAKData;
            close all;
            disp([rep, sub]);
        end

        %% Step 8: Concatenate Repetitions and Save Data
        for block = 1:3
            repN = R.import.repsel{subI};
            ftdata = appendFTData(ftdata_rep, repN, block);
            ftdata.hdr.label = ftdata.label;
            ftdata.fsample = ftdata_rep{R.import.repsel{subI}(1), block}.fsample;
            ftdata = makeContinuousTime(ftdata);

            % Add Electrode Locations
            load('acticap128_elec', 'elec');
            ftdata.elec = elec;

            % Save Data Locally
            ftdata = addHistoryField(ftdata, 'Import');
            if isfield(ftdata, 'cfg')
                ftdata = rmfield(ftdata, 'cfg');
            end
            saveExpData(R, sub{1}, blockpart{block}, [], 'rawdata', ftdata, 'raw');
        end

        %% Step 9: Save Unity Data
        UNData_cat = vertcat(UNData{:});
        saveExpData(R, sub{1}, 'Unity', [], '', UNData_cat, 'table');
        clear UNData UNData_cat;
    end
end
