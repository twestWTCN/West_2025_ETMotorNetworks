function enrichTrialInfo(R,subsel)
subI = 0;
for sub = subsel
    subI = subI + 1;

    % Correct Blocked Task Data
    corDat = [];
    corDat.sub =sub;
    corDat.block = 3; % task
    corDat.cond = 1:4; % condition
    corDat.part = 3; %2:5; % movement
    corDat.fileapend = "[R.epoch.names{part} '_trans_pp_arV3']";
    corDat.subfold = "['Condition' num2str(cond) '_epoched']";
    corDat.corrfx = @computeMoveTrialInfo;
    corDat.vargin = {R,sub,subI,0,'part'};

    correctData(R,corDat)

end