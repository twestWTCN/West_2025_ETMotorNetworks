function TAKData = importMotiveData(filename,fsamp,takconvpath,keycodepath,visSplit)
%% Convert the Motive TAK File to CSV
% You need to have the Motive converter installed:
%https://github.com/ha5dzs/optitrack-motive-file-converter
fileout = [filename(1:end-3) 'csv'];
system([takconvpath ' ' filename ' ' fileout]); % this runs the converter as a system cmd

%% Read in
% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 18);

% Specify range and delimiter
opts.DataLines = [8, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "TimeSeconds", "RotX", "RotY", "RotZ", "RotW", "X", "Y", "Z", "VarName10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"];
opts.SelectedVariableNames = ["TimeSeconds", "RotX", "RotY", "RotZ", "RotW", "X", "Y", "Z"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18"], "EmptyFieldRule", "auto");

% Import the data
TAKData = readtable(fileout, opts);
clear opts

%% Setup file

time_TAK = TAKData.TimeSeconds;
data_TAK = [TAKData.RotX TAKData.RotY TAKData.RotZ TAKData.RotW TAKData.X TAKData.Y TAKData.Z];

% Remove incomplete rows
data_TAK(isnan(time_TAK),:) = [];
time_TAK(isnan(time_TAK),:) = [];

time_TAK(any(isnan(data_TAK),2),:) = [];
data_TAK(any(isnan(data_TAK),2),:) = [];

% Resample to Target
tartime_TAK = time_TAK(1):1/fsamp:time_TAK(end);
% Ensure no duplicate points
[C,IA] = unique(time_TAK);
time_TAK = time_TAK(IA);
data_TAK = data_TAK(IA,:);

trials = interp1(time_TAK,data_TAK,tartime_TAK');
time_TAK = tartime_TAK;

%% Truncate
if visSplit
    figure
    plot(trials)
    xstart = ginput(2);
    hold on
    scatter(xstart(1,1),xstart(1,2),'g')
    hold on
    scatter(xstart(2,1),xstart(2,2),'r')
    
    trials = trials(xstart(1,1):xstart(2,1),:);
    time_TAK = time_TAK(xstart(1,1):xstart(2,1));

    delete(gcf)
end
%% Setup FT data structure
TAKData = [];
TAKData.fsample = fsamp;
TAKData.trial{1} = trials';
TAKData.time{1} = time_TAK;
TAKData.label = {'RotX','RotY','RotZ','RotW','xFinger1Pos','yFinger1Pos','zFinger1Pos'};
TAKData.chantype = {'MotiveRot','MotiveRot','MotiveRot','MotiveRot','MotivePos','MotivePos','MotivePos'};

%% Add the Screen Positions
XKeySaved = ReadCalibKeys([keycodepath 'XKey.txt']);
YKeySaved = ReadCalibKeys([keycodepath 'YKey.txt']);

if max(abs(XKeySaved))>1e3
    XKeySaved = ReadCalibKeysUK([keycodepath 'XKey.txt']);
    YKeySaved = ReadCalibKeysUK([keycodepath 'YKey.txt']);
end

trackingx = TAKData.trial{1}(strcmp('zFinger1Pos',TAKData.label),:)./1000;
trackingy = TAKData.trial{1}(strcmp('xFinger1Pos',TAKData.label),:)./1000;

xScreen = ((trackingx - XKeySaved(1)) * XKeySaved(2)) + XKeySaved(3);
yScreen = ((trackingy - YKeySaved(1)) * YKeySaved(2)) + YKeySaved(3);

TAKData.trial{1} = [TAKData.trial{1}; TAKData.trial{1}(end-2:end,:); xScreen; yScreen];
TAKData.label = {'RotX','RotY','RotZ','RotW','PosX','PosY','PosZ','xFinger1Pos','yFinger1Pos','zFinger1Pos','xScreen','yScreen'};
TAKData.chantype = {'MotiveRot','MotiveRot','MotiveRot','MotiveRot','MotivePos','MotivePos','MotivePos','MotivePos','MotivePos','MotivePos','UnityPos','UnityPos'};

save([fileparts(filename) '\TAKData_' filename(end-4) '.mat'],'TAKData')
