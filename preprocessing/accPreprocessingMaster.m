function data_out = accPreprocessingMaster(R, data_in)
    % ACCPREPROCESSINGMASTER: Preprocesses ACC/Motive data by applying low-pass filtering.
    % Parameters:
    %   R - Structure containing processing options (not used in the function currently)
    %   data_in - FieldTrip data structure containing ACC/Motive data
    % Returns:
    %   data_out - Processed ACC/Motive data

    % Initialize output as input
    data_out = data_in;

    % Select ACC/Motive channels
    accSel = strcmp(data_in.hdr.chantype, 'MotivePos') | strcmp(data_in.hdr.chantype, 'Accelerometer') | strcmp(data_in.hdr.chantype, 'UnityPos');

    %% Step 1: Low-Pass Filter (at 20 Hz)
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 20;
    cfg.lpfilttype = 'firws';
    cfg.lpinstabilityfix = 'reduce';  % Fix potential instability during filtering
    cfg.channel = data_in.label(accSel);
    data_in = ft_preprocessing(cfg, data_in);

    %% Step 2: Replace Original Data with Processed Data for ACC Channels
    for tr = 1:numel(data_out.trial)
        data_out.trial{tr}(accSel, :) = data_in.trial{tr};
    end
end