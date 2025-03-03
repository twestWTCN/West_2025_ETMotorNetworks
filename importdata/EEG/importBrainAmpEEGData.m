function EEGData = importBrainAmpEEGData(R,BAFN,visSplit)
if nargin<3
    visSplit = 2;
end

if strcmp(BAFN(end-2:end),'eeg')
    cfg = [];
    cfg.dataset = BAFN;
    EEGData = ft_preprocessing(cfg);
else
    load(BAFN,'trialData')
    EEGData = trialData;
end


% Resample to target
% % Memory inefficient:
% % cfg = [];
% % cfg.resamplefs = R.import.rs_fs;
% % EEGData = ft_resampledata(cfg,EEGData);

% Perform channel by channel resample (memory efficient)
datr = {};
for i=1:numel(EEGData.label)
    cfgp         = [];
    % instead of specifying channel names, you are allowed to use channel numbers
    cfgp.channel = i;
    datp         = ft_preprocessing(cfgp,EEGData);
    
    cfgr            = [];
    cfgr.resamplefs = R.import.rs_fs;
    datr{i}         = ft_resampledata(cfgr, datp);
%     clear datp
end

cfg = [];
EEGData = ft_appenddata(cfg, datr{:}); % this expands all cells into input variables
clear datr
% end
% Set data types
index = find(strncmp(EEGData.label,'Acc',3));
EEGData.chantype(index) = repmat({'Accelerometer'},numel(index),1);

index = find(strncmp(EEGData.label,'EMG',3));
EEGData.chantype(index) = repmat({'EMG'},numel(index),1);

index = find(strncmp(EEGData.label,'EKG',3));
EEGData.chantype(index) = repmat({'ECG'},numel(index),1);

index = find(strncmp(EEGData.label,'LabJack',3));
EEGData.chantype(index) = repmat({'Coder'},numel(index),1);

EEGData.chantype(1:128) = repmat({'EEG'},numel(index),1);

if visSplit
    [path,filen] = fileparts(BAFN);
    splitBA_EEG_byHand(EEGData,path,filen)
end


