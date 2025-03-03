function importExperimentalData_OPM(R)
    % IMPORTEXPERIMENTALDATA_OPM: Imports OPM, AV, and motion tracking data,
    % synchronizes, pads to same length, and saves the processed data.
    %   R - Structure containing paths and import options

    %% Initialize Variables
    close all;
    blockpart = R.import.blockpart;
    subI = 0;

    for sub = R.import.subsel
        subI = subI + 1;
        ftdata_rep = struct([]); % Struct to concatenate repetitions

        %% Iterate through Selected Repetitions
        for rep = R.import.repsel{subI}
            % Load AV Data (Audio-Visual)
            AVFN = [R.path.datapath '\' R.path.expname '\' R.import.subsel{subI} ...
                '\StimuliPCLocal\' sub{1} '_Rep_' num2str(rep) '\Data.txt'];
            try
                AVData = importEngAVData(AVFN, R.import.rs_fs);
            catch
                AVFN = [R.path.datapath '\' R.path.expname '\' R.import.subsel{subI} ...
                    '\StimuliPCLocal\' sub{1} '_Rep_' num2str(rep) '\Coder.txt'];
                AVData = importEngAVData_CODER(AVFN, R.import.rs_fs);
            end

            % Load Motion Tracking Data from MOTIVE
            TAKFN = [R.path.datapath '\' R.path.expname '\' R.import.subsel{subI} ...
                '\StimuliPCLocal\Motive\' sub{1} '_Rep_' num2str(rep) '.tak'];
            keycodepath = [R.path.datapath '\' R.path.expname '\' sub{1} '\StimuliPCLocal\' sub{1} '_Calibration\'];

            if ~exist([fileparts(TAKFN) '\TAKData_' num2str(rep) '.mat'], 'file')
                visSplit = 0;
                TAKData = importMotiveData(TAKFN, R.import.rs_fs, R.path.takconvpath, keycodepath, visSplit);
            else
                load([fileparts(TAKFN) '\TAKData_' num2str(rep) '.mat'], 'TAKData');
            end

            % Load OPM Data (MEG)
            opmChFN = [R.path.datapath '\' R.path.expname '\' R.import.subsel{subI} ...
                '\OPM\sub-' sub{1} '\ses-001\meg\sub-' sub{1} '_ses-001_task-' ...
                R.path.taskname{subI} '_run-00' num2str(rep) '_channels.tsv'];
            opmFN = [R.path.datapath '\' R.path.expname '\' R.import.subsel{subI} ...
                '\OPM\sub-' sub{1} '\ses-001\meg\sub-' sub{1} '_ses-001_task-' ...
                R.path.taskname{subI} '_run-00' num2str(rep) '_meg.json'];
            opmCSV = [R.path.datapath '\' R.path.expname '\' R.import.subsel{subI} '\HeadCast'];

            OPMData = importOPMData(R, opmFN, opmChFN, opmCSV, rep);

            %% Pad Data to Uniform Length
            datSize = [numel(AVData.time{1}), numel(OPMData.time{1}), numel(TAKData.time{1})];
            [~, dL] = max(datSize);
            timers = {AVData.time, OPMData.time, TAKData.time};

            % Pad AV, OPM, TAK data accordingly to match length
            if numel(AVData.time{1}) < datSize(dL)
                dSize = datSize(dL) - numel(AVData.time{1});
                AVData = padData(AVData, dSize, timers{dL}, 'randn');
            end
            if numel(OPMData.time{1}) < datSize(dL)
                dSize = datSize(dL) - numel(OPMData.time{1});
                OPMData = padData(OPMData, dSize, timers{dL}, 'randn');
            end
            if numel(TAKData.time{1}) < datSize(dL)
                dSize = datSize(dL) - numel(TAKData.time{1});
                TAKData = padData(TAKData, dSize, timers{dL}, 'randn');
            end

            %% Realign Data
            time_OPM = OPMData.time{1};
            trig_OPM = OPMData.trial{1}(strncmp(OPMData.label, 'Trigger_OPM', 3), :);
            time_AV = AVData.time{1};
            trig_AV = AVData.trial{1}(strncmp(AVData.label, 'Coder', 3), :);
            tarInd = [find(strncmp(OPMData.label, 'Trigger_OPM', 3)), find(strncmp(AVData.label, 'Coder', 3))];
            [OPMData.trial{1}, AVData.trial{1}] = realignData(time_OPM, trig_OPM, time_AV, trig_AV, R.import.rs_fs, OPMData.trial{1}, AVData.trial{1}, 1, tarInd);

            % Set OPM and AV data lengths to the minimum
            [OPMData, AVData] = setFTLengthToMin(OPMData, AVData);

            %% Concatenate OPM and AV Data
            chantypeConcat = [OPMData.chantype, AVData.chantype];
            OPMData = rmfield(OPMData, 'chantype');
            AVData = rmfield(AVData, 'chantype');
            OPMData = rmfield(OPMData, 'sampleinfo');
            AVData = rmfield(AVData, 'sampleinfo');

            cfg = [];
            ftdata = ft_appenddata(cfg, OPMData, AVData);
            ftdata.fsample = R.import.rs_fs;
            ftdata.hdr.chantype = chantypeConcat;
            ftdata.hdr.label = ftdata.label;

            %% Concatenate Repetitions and Save Data
            for block = 1:3
                cfg = [];
                cfg.trl = reshape(SMRtrialdef{block, 2}(:), [], 2);
                cfg.trl = [cfg.trl, zeros(size(cfg.trl, 1), 1)];
                ftdata_rep{rep, block} = ft_redefinetrial(cfg, ftdata);
            end

            clear ftdata_tmp AVData OPMData TAKData;
        end

        %% Save Data for All Blocks
        for block = 1:3
            repN = R.import.repsel{subI};
            ftdata = appendFTData(ftdata_rep, repN, block);
            ftdata.hdr.label = ftdata.label;
            ftdata.grad = gradStore;
            ftdata.fsample = ftdata_rep{R.import.repsel{subI}(1), block}.fsample;
            ftdata = makeContinuousTime(ftdata);
            ftdata = addHistoryField(ftdata, 'Import');
            if isfield(ftdata, 'cfg')
                ftdata = rmfield(ftdata, 'cfg');
            end
            saveExpData(R, sub{1}, blockpart{block}, [], 'rawdata', ftdata, 'raw');
        end
    end
end
