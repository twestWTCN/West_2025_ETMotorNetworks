function data_out = emgPreprocessingMaster(R, data_in)
    % EMGPREPROCESSINGMASTER: Preprocesses EMG data by applying high-pass filtering,
    % rectification, and low-pass filtering.
    % Parameters:
    %   R - Structure containing processing options (not used in the function currently)
    %   data_in - FieldTrip data structure containing EMG data
    % Returns:
    %   data_out - Processed EMG data

    % Initialize output as input
    data_out = data_in;

    % Select EMG channels
    emgSel = strncmp(data_in.label, 'EMG', 2);

    %% Step 1: High-Pass Filter (at 10 Hz)
    cfg = [];
    cfg.hpfilter = 'yes';
    cfg.hpfreq = 10;
    cfg.hpfilttype = 'firws';
    cfg.demean = 'yes';
    % Remove powerline interference (50 Hz, 100 Hz, 150 Hz)
    cfg.dftfreq = [50, 100, 150];
    cfg.dftreplace = 'neighbour';
    cfg.dftbandwidth = [1, 1, 1];
    cfg.dftneighbourwidth = [1, 1, 1];
    cfg.channel = data_in.label(emgSel);
    cfg.rectify = 'yes';  % Rectify the EMG signal
    data_in = ft_preprocessing(cfg, data_in);

    %% Step 2: Low-Pass Filter (at 48 Hz)
    cfg = [];
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 48;
    cfg.lpfilttype = 'firws';
    data_in = ft_preprocessing(cfg, data_in);

    %% Step 3: Replace Original Data with Processed Data for EMG Channels
    for tr = 1:numel(data_out.trial)
        data_out.trial{tr}(emgSel, :) = data_in.trial{tr};
    end
end
