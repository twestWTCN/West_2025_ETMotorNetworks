function AVData = importLocalDataGerman(filename,fsamp)
% Imports Unity data from Germany - 
pathstr = fileparts(filename);
AVData = readGermanAVData([pathstr '\RawData.mat']);

% Import the data
% AVData = readGermanDataFormat(filename);

time_AV = AVData.Timer;
data_AV = [AVData.Condition AVData.Coder AVData.xScreenPosition AVData.yScreenPosition AVData.trackingX AVData.trackingY AVData.trackingZ];

% Remove incomplete rows
data_AV(isnan(time_AV),:) = [];
time_AV(isnan(time_AV),:) = [];

time_AV(any(isnan(data_AV),2),:) = [];
data_AV(any(isnan(data_AV),2),:) = [];

% Resample to Target
tartime_AV = time_AV(1):1./fsamp:time_AV(end);
% Ensure no duplicate points
[C,IA] = unique(time_AV);
time_AV = time_AV(IA);

data_AV = data_AV(IA,:);

trials = interp1(time_AV,data_AV,tartime_AV');
time_AV = tartime_AV;

%% Setup FT data structure
AVData = [];
AVData.fsample = fsamp;
AVData.trial{1} = trials';
AVData.time{1} = time_AV;
AVData.label = {'Condition','Coder','xScreen','yScreen','xFinger1Pos','yFinger1Pos','zFinger1Pos'};
AVData.chantype = {'Info','Trigger','UnityPos','UnityPos','MotivePos','MotivePos','MotivePos'};
% Save to file
save([fileparts(filename) '\AVData.mat'],'AVData')

end