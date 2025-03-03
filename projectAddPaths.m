function R = projectAddPaths(R)
switch getenv('computername')
    case 'IC-H5TFLZ3'
        gitpath = 'Z:\TimWest\GITHUB';
        R.path.spmpath = '\\192.168.0.10\cagnan_lab\TimWest\GITHUB\spm12';
        R.path.datapath = '\\192.168.0.10\cagnan_lab\TimWest\DATA\Data';
        R.path.spm12path = '\\192.168.0.10\cagnan_lab\TimWest\GITHUB\spm12';
        R.path.spm8path = '\\192.168.0.10\cagnan_lab\TimWest\GITHUB\spm8';
        R.path.ftpath = '\\192.168.0.10\cagnan_lab\TimWest\GITHUB\fieldtrip\';
        R.path.takconvpath = 'C:\TAKconverter\converter.exe';
        
        %% Add your computer name as
        % case '#COMPUTERNAME#'
        % R.path.spmpath = ''; % SPM default path
        % R.path.datapath = '';% Data path
        % R.path.spm12path = ''; % SPM12  path
        % R.path.spm8path = ''; % SPM8  path
        % R.path.ftpath = ''; % FieldTrip Path
        % R.path.takconvpath = ''; % TAK converter path
end

pathCell = regexp(path, pathsep, 'split'); onPath = any(strcmpi(R.path.spm8path, pathCell)); %import only
if ~onPath; addpath(R.path.spm12path); spm eeg; close all; end; %import only

rmpath(genpath(fullfile(R.path.spm12path,'external','fieldtrip'))); %import only
addpath(R.path.ftpath); ft_defaults();

