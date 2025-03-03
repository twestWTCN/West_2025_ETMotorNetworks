function preprocessContinuousData(R)
    % PREPROCESSCONTINUOUSDATA: Preprocesses continuous EEG data by importing,
    % cleaning, and performing ICA-based artifact removal.
    % Parameters:
    %   R - Structure containing paths and preprocessing options

    % Load default electrode configuration
    elec_default = ft_read_sens([R.path.ftpath '\template\electrode\standard_1005.elc']);
    elec_default = ft_convert_units(elec_default, 'mm');

    %% Step 1: Main Loop to Process Each Subject
    blockpart = R.import.blockpart;  % Specify parts to import
    subI = 0;

    for sub = R.preprocess.subsel  % Processing selected subjects
        subI = subI + 1;

        % Load experimental data for each block part (Rest, Posture, Task)
        ftdata{1} = loadExpData(R, sub{1}, blockpart{1}, [], 'rawdata', 'raw');
        samptype = zeros(1, numel(ftdata{1}.trial));

        ftdata{2} = loadExpData(R, sub{1}, blockpart{2}, [], 'rawdata', 'raw');
        samptype = [samptype, ones(1, numel(ftdata{2}.trial))];

        ftdata{3} = loadExpData(R, sub{1}, blockpart{3}, [], 'rawdata', 'raw');
        samptype = [samptype, repmat(2, 1, numel(ftdata{3}.trial))];

        % Concatenate all blocks into a single dataset
        ftdata_cat = ft_appenddata([], ftdata{:});
        ftdata_cat.hdr = ftdata{1}.hdr;  % Assign header information

        % Assign electrode positions if the site is DL
        if strcmp(R.import.site, 'DL')
            ftdata_cat.elec = elec_default;
        end
        clear ftdata;

        %% Step 2: Basic Preprocessing (including ICA)
        [ftdata_cat, ~, badTr] = basicPreProcessing(R, ftdata_cat);
        samptype = samptype(~badTr);  % Correct sample type for removed trials

        %% Step 3: Split Data Back into Blocks and Save
        for block = 1:numel(blockpart)
            cfg = [];
            cfg.trials = find(samptype == block - 1);  % Select trials for the specific block
            ftdata = ft_selectdata(cfg, ftdata_cat);

            % Check header consistency
            checkHeader(ftdata);

            % Add preprocessing history and save the data
            ftdata = addHistoryField(ftdata, 'Preprocessed');
            fileappend = 'pp';
            saveExpData(R, sub{1}, blockpart{block}, [], fileappend, ftdata, 'pp');
        end

        clear ftdata;
    end
end
