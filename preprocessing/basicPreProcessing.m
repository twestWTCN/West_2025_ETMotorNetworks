function [ftdata_pp,R,badTr] = basicPreProcessing(R,ftdata,code)
ftdata_pp = ftdata;
%% Replace NaNs with zeros
for tr = 1:numel(ftdata_pp.trial)
    X = ftdata_pp.trial{tr};
    X(isnan(X)) = 0;
    ftdata_pp.trial{tr} = X;
end

% Now go through modalities
if any(strcmp(ftdata_pp.hdr.chantype,'MotivePos')) || any(strcmp(ftdata_pp.hdr.chantype,'UnityPos')) || any(strcmp(ftdata_pp.hdr.chantype,'Accelerometer'))
    % ACC/MOT
    ftdata_pp = accPreprocessingMaster(R,ftdata_pp);
else
    warning('No ACCs/Mots found, so not preprocessing!')
end

if any(strncmp(ftdata_pp.label,'emg',3)) || any(strncmp(ftdata_pp.label,'EMG',3))
    % EMG
    ftdata_pp = emgPreprocessingMaster(R,ftdata_pp);
else
    warning('No EMGs found, so not preprocessing!')
end
% OPM Preprocessing
if any(strncmp(ftdata_pp.label,'G',1))
    [ftdata_pp,badTr] = opmPreprocessingMaster(R,ftdata_pp,1);
else
    warning('No OPMs found, so not preprocessing!')
end

% EEG Preprocessing
if any(strncmp(ftdata_pp.hdr.chantype,'EEG',3) | strncmp(ftdata_pp.hdr.chantype,'eeg',3))
    [ftdata_pp,badTr] = eegPreprocessingMaster(R,ftdata_pp,1);
else
    warning('No EEGs found, so not preprocessing!')
end
