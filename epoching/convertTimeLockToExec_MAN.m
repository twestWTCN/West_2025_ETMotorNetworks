function [exec_trl, hold_trl, react_time, testdata] = convertTimeLockToExec_MAN(ftdata_tmp, trl, tranSamp, plotop, accBadFlag, searchwin)
% CONVERTTIMELOCKTOEXEC_MAN converts stimulus-locked trials to movement onset markers.
% This function extracts data, applies filters, and allows manual annotation of movement onset using user inputs.

close all;
fsamp = ftdata_tmp.fsample;
coder = ftdata_tmp.trial{1}(strncmp(ftdata_tmp.label, 'Coder', 4), :);

% Determine data channels to use (MotivePos or Accelerometer based on accBadFlag)
if accBadFlag == 1
    % Use MotivePos channels only
    motInd = find(strncmp(ftdata_tmp.hdr.chantype, 'MotivePos', 4));
else
    % Use Accelerometer channels along with MotivePos
    accInd = find(strncmp(ftdata_tmp.hdr.chantype, 'Accelerometer', 4));
    motInd = find(strncmp(ftdata_tmp.hdr.chantype, 'MotivePos', 4) & strncmp(ftdata_tmp.label, 'Pos', 3));
    motInd = union(motInd, find(strncmp(ftdata_tmp.label, 'xFinger1Pos', 4) | strncmp(ftdata_tmp.label, 'yFinger1Pos', 4) | strncmp(ftdata_tmp.label, 'zFinger1Pos', 4)));
    motInd(motInd == 0) = [];
end

% Extract accelerometer data if not using bad flag
if ~accBadFlag
    cfg = [];
    cfg.channel = ftdata_tmp.label(accInd);
    ftdataA = ft_selectdata(cfg, ftdata_tmp); % Select accelerometer data
end

% Extract motive data
cfg = [];
cfg.channel = ftdata_tmp.label(motInd);
ftdata_tmp = ft_selectdata(cfg, ftdata_tmp); % Select motive data
XD = ftdata_tmp.trial{1};

% Apply low-pass filter to reduce tremor
fc = 1;
[b, a] = butter(4, fc / (fsamp / 2));
XD = filtfilt(b, a, XD')';
XD = normalize(XD')';

% Use accelerometer data for further refinement if available
if ~accBadFlag
    accdata = ftdataA.trial{1};
    [b, a] = butter(2, fc / (fsamp / 2));
    accdata = filtfilt(b, a, accdata')';
    accdata = normalize(accdata')';
    clear ftdataA;
end

figure(10);
set(gcf, 'Position', [859, 90, 755, 798]);
movtime = []; react_time = []; markbank = []; ttestdata = [];
tri = 0;

% Loop through each trial and manually annotate movement
while tri < size(trl, 1)
    tri = tri + 1;
    if (trl(tri, 1) > 0) && (trl(tri, 2) <= size(ftdata_tmp.time{1}, 2))
        % Extract movement data from trial
        mdata = XD(:, trl(tri, 1):trl(tri, 2));
        mdata = mdata - median(mdata, 2);
        
        if ~accBadFlag
            adata = accdata(:, trl(tri, 1):trl(tri, 2));
            adata = adata - median(adata, 2);
        else
            adata = [];
        end
        
        ttestdata(:, :, tri) = [mdata; adata];
        tvec = linspace(searchwin(1) / fsamp, searchwin(2) / fsamp, size(mdata, 2));
        
        flag = 1;
        while flag
            subplot(2, 1, 1);
            a(1) = gca;
            plot(tvec, mdata(1:3, :), 'r');
            hold on;
            plot(tvec, adata, 'b');
            XR = [mdata; adata];
            XR1 = [min(XR(:)), max(XR(:))];
            title(num2str(tri / size(trl, 1)));
            
            subplot(2, 1, 2);
            a(2) = gca;
            plot(tvec(2:end), diff(mdata(1:3, :), [], 2), 'r');
            hold on;
            plot(tvec(2:end), diff(adata, [], 2), 'b');
            XR = [diff(mdata, [], 2); diff(adata, [], 2)];
            XR2 = [min(XR(:)), max(XR(:))];
            
            % User manual input for movement onset markers
            trtyp = {'r--', 'b--'};
            movlocal = [];
            for i = 1:2
                [x, ~] = ginput(1);
                subplot(2, 1, 1);
                plot([x(1), x(1)], XR1, trtyp{i});
                subplot(2, 1, 2);
                plot([x(1), x(1)], XR2, trtyp{i});
                [~, movlocal(i)] = min((tvec - x).^2); % Find the sample closest to the time stamp
            end
            
            disp('If happy press Space, if bad press r; else press another key to try again');
            w = waitforbuttonpress;
            kin = 'a';
            if w
                kin = get(gcf, 'CurrentCharacter');
            end
            if kin == '1' || kin == '2'
                flag = 0;
                badtr(i) = 0;
                movtime{1}(tri) = trl(tri, 1) + movlocal(1); % In samples
                movtime{2}(tri) = trl(tri, 1) + movlocal(2);
                
                react_time{1}(tri) = movlocal(1) / fsamp;
                react_time{2}(tri) = movlocal(2) / fsamp;
                rating(tri) = str2double(kin);
                markbank(:, tri) = movlocal;
                
            elseif kin == '3'
                movtime{1}(tri) = nan;
                movtime{2}(tri) = nan;
                
                react_time{1}(tri) = nan;
                react_time{2}(tri) = nan;
                flag = 0;
                badtr(i) = 1;
                rating(tri) = 3;
                markbank(:, tri) = [-1, -1];
            elseif kin == ',' || isempty(kin) % Go back
                tri = tri - 2;
                flag = 0;
            end
        end
        clf;
    else
        % Handle trials longer than available data
        warning('Trial definition is longer than recording, replacing with nans');
        movtime{1}(tri) = nan;
        movtime{2}(tri) = nan;
        react_time{1}(tri) = nan;
        react_time{2}(tri) = nan;
        markbank(:, tri) = [-1, -1];
        rating(tri) = nan;
    end
end

% Store the raw data as test
if size(ttestdata, 3) ~= size(markbank, 2)
    warning('Test data does not match the markbank');
end
if size(ttestdata, 3) ~= size(react_time{1}, 2)
    warning('Mismatch between test data and reaction time');
end

testdata.data = ttestdata;
testdata.trialmarker = markbank;
testdata.rating = rating;
testdata.reactbank = react_time;

% Shave off bad trials
trl(isnan(movtime{1}), :) = [];
react_time{1}(isnan(movtime{1})) = [];
react_time{2}(isnan(movtime{1})) = [];
rating(isnan(movtime{1})) = [];
movtime{1}(isnan(movtime{2})) = [];
movtime{2}(isnan(movtime{2})) = [];

% Redefine exec and hold trials
exec_trl = [movtime{1} + tranSamp(1); movtime{1} + tranSamp(2); repmat(tranSamp(1), 1, size(movtime{1}, 2)); trl(:, 4:end)'; react_time{1}; react_time{2}; rating]';
hold_trl = [movtime{2} + tranSamp(1); movtime{2} + tranSamp(2); repmat(tranSamp(1), 1, size(movtime{2}, 2)); trl(:, 4:end)'; react_time{1}; react_time{2}; rating]';

newtrl = exec_trl;

% Plot results if plotop is enabled
if plotop == 1
    ax(1) = subplot(2, 1, 1);
    plot(ftdata_tmp.time{1}, ftdata_tmp.trial{1}(: , :)', 'LineWidth', 2);
    ax(2) = subplot(2, 1, 2);
    X = ftdata_tmp.trial{1};
    mX = [min(X(:)), max(X(:))];
    a = plot(ftdata_tmp.time{1}, X', 'LineWidth', 2);
    yyaxis right;
    plot(ftdata_tmp.time{1}, coder', 'LineWidth', 1);
    yyaxis left;
    for tri = 1:size(newtrl, 1)
        try
            hold on;
            plot([ftdata_tmp.time{1}(newtrl(tri, 1)), ftdata_tmp.time{1}(newtrl(tri, 2))], [mX(1), mX(1)], 'k--');
            plot([ftdata_tmp.time{1}(newtrl(tri, 1)), ftdata_tmp.time{1}(newtrl(tri, 2))], [mX(2), mX(2)], 'k--');
            plot([ftdata_tmp.time{1}(newtrl(tri, 1)), ftdata_tmp.time{1}(newtrl(tri, 1))], mX, 'k--');
            plot([ftdata_tmp.time{1}(newtrl(tri, 1) - tranSamp(1)), ftdata_tmp.time{1}(newtrl(tri, 1) - tranSamp(1))], mX, 'k-');
            plot([ftdata_tmp.time{1}(newtrl(tri, 2)), ftdata_tmp.time{1}(newtrl(tri, 2))], mX, 'k--');
        catch
            disp('Trial definition probably exceeds length of data!');
        end
    end
    legend(a, ftdata_tmp.label);
    xlabel('Time (s)'); ylabel('Movement (mm)');
    linkaxes(ax, 'x');
end

end
