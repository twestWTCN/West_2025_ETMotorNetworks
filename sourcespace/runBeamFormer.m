function [source, frqtlockdataCat] = runBeamFormer(ftdata, headmodel, leadfield, rawtrialFlag, srcFilter, type, band, modality, refchan, frqtlockdataCat, bstrpflag)
% Run beamformer analysis
% Input: 
%   ftdata - FieldTrip data structure
%   headmodel - head model for beamforming
%   leadfield - lead field model
%   rawtrialFlag - flag for raw trial analysis
%   srcFilter - source filter
%   type - type of beamformer (DICS, LCMV, etc.)
%   band - frequency band for analysis
%   modality - EEG or OPM
%   refchan - reference channel (for DICSREF)
%   frqtlockdataCat - optional frequency-locked data concatenation
%   bstrpflag - bootstrap flag

if nargin < 10
    frqtlockdataCat = [];
end
if nargin < 11
    bstrpflag = 0;
end

%% Estimate the truncation point for the kappa
cfg = [];
cfg.covariance = 'yes';
switch modality
    case 'EEG'
        cfg.channel = 'eeg';
    case 'OPM'
        cfg.channel = ftdata.label(strncmp(ftdata.label, 'G', 1));
end
cfg.covariancewindow = 'all';
tlock = ft_timelockanalysis(cfg, ftdata);
S = svd(tlock.cov);

% Estimate kappa value based on modality
switch modality
    case 'OPM'
        [~, s, ~] = svd(tlock.cov);
        d = -diff(log10(diag(s)));
        d = d ./ std(d);
        kappaEst = find(d(10:end) > 2, 1, 'first') + 10;
    case 'EEG'
        [~, s, ~] = svd(tlock.cov);
        d = -diff(log10(diag(s)));
        d = d ./ std(d);
        kappaEst = find(d(10:end) > 2, 1, 'first') + 10;
end
kappaEst = kappaEst - 1;

if isempty(kappaEst)
    kappaEst = length(S);
end

h = figure;
plot(1:length(S), log10(S));
hold on;
plot([kappaEst kappaEst], [min(log10(S)) max(log10(S))], 'r');
box off;
shg;
pause(2); close(h);

%% Compute the beamformers
switch type
    case {'DICS', 'ELO', 'DICSREF', 'DICSREFI'}
        if isempty(frqtlockdataCat)
            trl = 0; seg = 0;
            while trl < numel(ftdata.trial)
                seg = seg + 1; trl = trl + 1;
                % Break up into manageable batches
                cfg = [];
                if trl + 50 < numel(ftdata.trial)
                    cfg.trials = trl:trl + 50;
                else
                    cfg.trials = trl:numel(ftdata.trial);
                end
                ftdata_brk = ft_selectdata(cfg, ftdata);
                trl = trl + 50;
                
                % Prepare Data
                cfg = [];
                if strcmp(type, 'DICS')
                    switch modality
                        case 'EEG'
                            cfg.channel = 'eeg';
                        case 'OPM'
                            cfg.channel = ftdata.label(strncmp(ftdata.label, 'G', 1));
                    end
                elseif strcmp(type, 'DICSREF') || strcmp(type, 'DICSREFI')
                    switch modality
                        case 'EEG'
                            cfg.channel = ['eeg' refchan];
                        case 'OPM'
                            cfg.channel = [ftdata.label(strncmp(ftdata.label, 'G', 1)); refchan];
                    end
                end
                
                cfg.method = 'mtmfft';
                cfg.pad = 'nextpow2';
                cfg.output = 'powandcsd';
                cfg.tapsmofrq = 3;
                cfg.foi = band(1):0.5:band(2);
                cfg.keeptrials = 'yes';
                freqdata = ft_freqanalysis(cfg, ftdata_brk);
                csd = freqdata.crsspctrm;
                
                if seg == 1
                    freqdataCat = freqdata;
                    csdCat = csd;
                else
                    cfg = [];
                    cfg.keepsampleinfo = 'no';
                    cfg.appenddim = 'rpt';
                    freqdataCat = ft_appendfreq(cfg, freqdataCat, freqdata);
                    csdCat(end+1:end+size(csd, 1), :, :) = csd;
                end
            end
            freqdataCat.labelcmb = freqdata.labelcmb;
            freqdataCat.crsspctrm = csdCat;
            delete csdCat freqdata;
        else
            freqdataCat = frqtlockdataCat;
        end
        
        cfg = [];
        if isempty(srcFilter)
            cfg.channel = intersect(leadfield.label, freqdataCat.label);
        else
            cfg.channel = intersect(intersect(leadfield.label, freqdataCat.label), srcFilter.label);
        end
        leadfield = ft_selectdata(cfg, leadfield);
        
        if ~isempty(srcFilter)
            leadfield.label = srcFilter.label;
        end
        
        % DICS Beamformer
        cfg = [];
        cfg.method = 'dics';
        cfg.frequency = median(band);
        cfg.sourcemodel = leadfield;
        cfg.headmodel = headmodel;
        cfg.dics.projectnoise = 'yes';
        cfg.dics.lambda = '5%';
        cfg.dics.kappa = kappaEst;
        if strcmp(type, 'DICSREFI')
            cfg.dics.powmethod = 'imag';
            warning('You are using imaginary coherence!');
        end
        if isempty(srcFilter)
            cfg.dics.keepfilter = 'yes';
        else
            cfg.sourcemodel.filter = srcFilter.filter;
        end
        cfg.keeptrials = 'yes';
        cfg.rawtrial = rawtrialFlag;
        if strcmp(type, 'DICSREF') || strcmp(type, 'DICSREFI')
            cfg.refchan = refchan;
        end
        if bstrpflag
            cfg.bootstrap = 'yes';
            cfg.numbootstrap = 200;
        else
            cfg.bootstrap = 'no';
        end
        source = ft_sourceanalysis(cfg, freqdataCat);
        if strcmp(type, 'DICS')
            if isempty(srcFilter)
                source.avg.pow = source.avg.pow ./ source.avg.noise;
            elseif strcmp(rawtrialFlag, 'yes')
                for i = 1:numel(source.trial)
                    source.trial(i).pow = source.trial(i).pow ./ source.trial(i).noise;
                end
            else
                source.avg.pow = source.avg.pow ./ source.avg.noise;
            end
        end
        frqtlockdataCat = freqdataCat;
        
    case 'LCMV'
        if isempty(frqtlockdataCat)
            trl = 0; seg = 0;
            while trl < numel(ftdata.trial)
                seg = seg + 1; trl = trl + 1;
                % Break up into manageable batches
                cfg = [];
                if trl + 25 < numel(ftdata.trial)
                    cfg.trials = trl:trl + 25;
                else
                    cfg.trials = trl:numel(ftdata.trial);
                end
                ftdata_brk = ft_selectdata(cfg, ftdata);
                trl = trl + 50;
                
                % Prepare Data
                cfg = [];
                cfg.covariance = 'yes';
                switch modality
                    case 'EEG'
                        cfg.channel = 'eeg';
                    case 'OPM'
                        cfg.channel = ftdata.label(strncmp(ftdata.label, 'G', 1));
                end
                cfg.keeptrials = 'yes';
                cfg.covariancewindow = 'all';
                tlock = ft_timelockanalysis(cfg, ftdata);
                
                if seg == 1
                    tlockCat = tlock;
                else
                    cfg = [];
                    cfg.keepsampleinfo = 'no';
                    cfg.appenddim = 'rpt';
                    tlockCat = ft_appendtimelock(cfg, tlockCat, tlock);
                end
            end
        else
            tlockCat = frqtlockdataCat;
        end
        
        cfg = [];
        if isempty(srcFilter)
            cfg.channel = intersect(leadfield.label, tlockCat.label);
        else
            cfg.channel = intersect(intersect(leadfield.label, tlockCat.label), srcFilter.label);
        end
        leadfield = ft_selectdata(cfg, leadfield);
        
        if ~isempty(srcFilter)
            leadfield.label = srcFilter.label;
        end        
        
        % LCMV Beamformer
        cfg = [];
        cfg.method = 'lcmv';
        cfg.sourcemodel = leadfield;
        cfg.headmodel = headmodel;
        cfg.lcmv.projectnoise = 'yes';
        cfg.lcmv.lambda = '1%';
        cfg.lcmv.kappa = kappaEst;
        cfg.lcmv.weightnorm = 'unitnoisegain';
        if isempty(srcFilter)
            cfg.lcmv.keepfilter = 'yes';
        else
            cfg.lcmv.filter = srcFilter.filter;
        end
        cfg.keeptrials = 'yes';
        cfg.rawtrial = rawtrialFlag;
        source = ft_sourceanalysis(cfg, tlockCat);
        
        if isempty(srcFilter)
            source.avg.pow = source.avg.pow ./ source.avg.noise;
        elseif strcmp(rawtrialFlag, 'yes')
            for i = 1:numel(source.trial)
                source.trial(i).pow = source.trial(i).pow ./ source.trial(i).noise;
            end
        else
            source.avg.pow = source.avg.pow ./ source.avg.noise;
        end
        frqtlockdataCat = tlockCat;
end

% Clear cfg field to save memory
source = rmfield(source, 'cfg');

end
