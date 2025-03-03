function checkEpochAcc(ftdata_trans,accBadFlag)
if accBadFlag == 1
    accInd = find(strncmp(ftdata_trans.hdr.chantype,'MotivePos',4));
else % use the Acceler position Only
%     accInd = find(strncmp(ftdata_trans.hdr.chantype,'Accelerometer',4)| strncmp(ftdata_trans.hdr.chantype,'MotivePos',4) );
%     accInd = find(strncmp(ftdata_trans.hdr.chantype,'MotivePos',4) );
    accInd = find(strncmp(ftdata_trans.hdr.chantype,'Accelerometer',4) );
end

cfg = [];
cfg.channel = ftdata_trans.label(accInd);
ftdata_trans = ft_selectdata(cfg,ftdata_trans);

% cfg = [];
% cfg.lpfreq = 'yes';
% cfg.bsfreq = [.5 10];
% ftdata_trans = ft_preprocessing(cfg,ftdata_trans);

clf
for i = 1:numel(ftdata_trans.trial)
    subplot(2,1,1)
    XD = ftdata_trans.trial{i};
    fc = .5;
    [b,a] = butter(4,fc/(512/2));
    XD = filtfilt(b,a,XD')';
    XD  = XD-mean(XD(1:256));
    plot(ftdata_trans.time{i},XD)
    hold on
    subplot(2,1,2)
    mdata = ftdata_trans.trial{i}(1:3,:);
    [dum,pmdat] = pca(abs(mdata)); % get main component
    mdata = pmdat.M(:,1)';
    %     mdata = (mdata-mean(mdata))./std(mdata);
    mdata  = mdata-mean(mdata(1:256));
%     mdata= mdata.^2;
%     mdata = mdata-median(mdata(:));
    plot(ftdata_trans.time{i},mdata)
    hold on
end