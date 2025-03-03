function OPMData = importOPMData(R, opmFN, opmChFN, opmCSV, repN)
    % IMPORTOPMDATA: Imports and processes OPM (Optically Pumped Magnetometers) data.
    % Converts the imported data to FieldTrip format and prepares it for further analysis.
    % Parameters:
    %   R - Structure containing paths and options
    %   opmFN - Filename for OPM MEG data
    %   opmChFN - Filename for OPM channel lists
    %   opmCSV - Path to CSV containing headcast positions
    %   repN - Repetition number for data loading

    %% Step 1: Generate Paths and Prepare Data
    path_to_data = [fileparts(opmFN) '\'];
    tsvpath = [opmFN(1:end-8) 'channels.tsv'];

    % Generate positions for hybrid headcasts
    try
        info2pos_TW2('BOTH', tsvpath, opmCSV);
    catch
        warning('Assuming that this headcast is a gen2/3 hybrid');
        info2pos_g3_g2_TW(tsvpath, opmCSV);
    end

    %% Step 2: Load OPM Data
    if exist([opmFN(1:end-4) 'dat'], 'file') > 0
        % If the data file already exists, load it
        D = spm_eeg_load([opmFN(1:end-4) 'mat']);
    else
        % Create a new data file if it doesn't exist
        S = [];
        S.data = [opmFN(1:end-4) 'bin'];
        S.channels = opmChFN;
        S.meg = opmFN;
        S.positions = [tsvpath(1:end-12) 'positions.tsv'];
        S.sMRI = [opmCSV '\mri.nii'];
        D = spm_opm_create(S);
    end

    %% Step 3: Convert Data to FieldTrip Format
    if exist(['OPMData_tmp_rep' num2str(repN) '.mat'], 'file') > 0
        load(['OPMData_tmp_rep' num2str(repN) '.mat'], 'OPMData');
    else
        OPMData = spm2fieldtrip(D);

        % Memory-efficient resampling of the data
        NSplit = 5;
        OPMData = memResample(OPMData, R.import.rs_fs, NSplit);
        save(['OPMData_tmp_rep' num2str(repN) '.mat'], 'OPMData');
    end

    %% Step 4: Remove Non-Gradient Channels
    fullChan = OPMData.label;
    nonGrad = setdiff(OPMData.label, OPMData.grad.label);
    for i = 1:numel(nonGrad)
        if ~strcmp(nonGrad{i}(1:2), 'NI')
            fullChan(strcmp(fullChan, nonGrad{i})) = [];
        end
    end

    cfg = [];
    cfg.channel = fullChan;
    OPMData = ft_selectdata(cfg, OPMData);

    % Remove unused gradient labels
    [~, IA] = setdiff(OPMData.grad.label, OPMData.label);
    OPMData.grad.label(IA) = [];
    OPMData.grad.chanori(IA, :) = [];
    OPMData.grad.chanpos(IA, :) = [];
    OPMData.grad.chanunit(IA, :) = [];
    OPMData.grad.coilpos(IA, :) = [];
    OPMData.grad.coilori(IA, :) = [];
    OPMData.grad.tra(IA, :) = [];
    OPMData.grad.tra(:, IA) = [];

    %% Step 5: Rename and Sort Channels
    % Rename Trigger channel
    OPMData.label{strcmp(OPMData.label, R.opmsense.trigchan{1})} = 'Trigger_OPM';

    % Rename Accelerometer channels
    accLab = {'AccX', 'AccY', 'AccZ'};
    for i = 1:3
        OPMData.label{strcmp(OPMData.label, R.opmsense.ACCchan{i})} = accLab{i};
    end

    %% Step 6: Define Channel Types
    OPMData.chantype = cell(1, size(OPMData.label, 1));
    index = ~cellfun('isempty', regexp(OPMData.label, 'Y$'));
    OPMData.chantype(index) = repmat({'radial'}, sum(index), 1);

    index = ~cellfun('isempty', regexp(OPMData.label, 'Z$'));
    OPMData.chantype(index) = repmat({'tangential'}, sum(index), 1);

    index = ~cellfun('isempty', regexp(OPMData.label, 'X$'));
    OPMData.chantype(index) = repmat({'axial'}, sum(index), 1);

    index = find(strncmp(OPMData.label, 'Sync', 4));
    OPMData.chantype(index) = repmat({'Noise'}, numel(index), 1);

    index = find(strncmp(OPMData.label, 'Flux', 4));
    OPMData.chantype(index) = repmat({'Flux'}, numel(index), 1);

    index = find(strncmp(OPMData.label, 'Trigger_OPM', 4));
    OPMData.chantype(index) = repmat({'Trigger'}, numel(index), 1);

    index = find(strncmp(OPMData.label, 'Acc', 3));
    OPMData.chantype(index) = repmat({'Accelerometer'}, numel(index), 1);

    %% Step 7: Remove Extra Trigger Channels
    index = strncmp(OPMData.label, 'NI-TRIG', 7);
    OPMData.chantype(index) = [];
    cfg = [];
    cfg.channel = OPMData.label(~index);
    OPMData = ft_selectdata(cfg, OPMData);
end
