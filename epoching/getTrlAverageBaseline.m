function baseMed = getTrlAverageBaseline(ftdata_tmp,trl,accBadFlag)
% Backcompatible with bad flag
motInd = find(strncmp(ftdata_tmp.hdr.chantype,'MotivePos',4));

cfg = [];
cfg.channel = ftdata_tmp.label(motInd);
ftdata_tmp = ft_selectdata(cfg,ftdata_tmp);
XD = ftdata_tmp.trial{1};

fs = ftdata_tmp.fsample;
fc = 2;
[b,a] = butter(4,fc/(fs/2));
XD = filtfilt(b,a,XD')';

% Convert from stimulus locked onset, to movement onset using window
baseMed = [];
for tri = 1:size(trl,2)
    baseMed(:,tri) = median(XD(:,trl(1,tri):trl(2,tri)),2);
end
baseMed = median(baseMed,2);