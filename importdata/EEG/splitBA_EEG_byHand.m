function splitBA_EEG_byHand(EEGData,path,filen)
tvec = EEGData.time{1};
coder = EEGData.trial{1}(strcmp(EEGData.label,'LabJack'),:);
plot(coder)

nEpoch = input('How many epochs?');

disp('Now start points')

xstart = ginput(nEpoch);
hold on
scatter(xstart(:,1),xstart(:,2),'g')

disp('Now end points')
xend = ginput(nEpoch);
hold on
scatter(xend(:,1),xend(:,2),'r')

for i = 1:nEpoch
    mark = fix(xstart(i,1));
    backtrack = fix(40.*EEGData.fsample);
    
    if (mark-backtrack)<1
        backtrack = mark-1;
    end
    XC = diff(coder(1,mark-backtrack:mark));
    
    if isempty(XC)
        XC = mark;
    else
        realStart = mark - (numel(XC)-find(abs(XC)>500,1,'last'));
    end
    
    if i ~= nEpoch
        
        mark = fix(xend(i,1));
        backtrack = fix(40.*EEGData.fsample);
        XC = diff(coder(1,mark-backtrack:mark));
        
        if isempty(XC)
            XC = mark;
        else
            realEnd = mark - (numel(XC)-find(abs(XC)>100,1,'last'));
        end
    else
        realEnd = size(tvec,2);
    end
    plot([realStart realStart],[min(coder) max(coder)],'g--')
    plot([realEnd realEnd],[min(coder) max(coder)],'r--')
    
    cfg = [];
    cfg.toilim = [tvec(realStart) tvec(realEnd)];
    trialData = ft_redefinetrial(cfg,EEGData);
    save([path '\' filen '_Rep' num2str(i)],'trialData');
    
end

