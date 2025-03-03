function [decoder,SMRtrialdef] = codetimingsEEG(CVM,ampCoder,trialdetails,splitType,blockSize,minbreak)

% ---- ampCoder = Coder coming from the amplifier of a movement/calibration trial (port #14)
% ---- CVM      = table extracted from trigger calibration measurement

% Try to decode using simple absolute differences
ampCoder(isnan(ampCoder)) = 0;
ampDeCoder = abs(ampCoder-CVM); % absolute Difference
[dError,decoder] = min(ampDeCoder,[],1);
decoder(dError>0.01) = nan;
% Now make a trialdef variable that you can use to segment the data in
% fieldtrip
for trialstage = 1:numel(CVM) % loop through number of trial stages
    X = find(decoder==trialstage); % find timings at which trial stage exists
    blocks = SplitVec(X,'consecutive');
    
    % Now loop through number of instances/stage repetitions SR
    nmx = [];
    for SR = 1:numel(blocks)
        nmx(SR) = numel(blocks{SR});
        trialdef{trialstage}(:,SR) = [blocks{SR}(1) blocks{SR}(end)]; % i.e. take first and last samples of each block
        % Make sure begin samples are not allocated to a trial stage:
    end
    
    %% Remove small breaks (minor collisions)
    trBreak = [inf trialdef{trialstage}(2,2:end)-trialdef{trialstage}(2,1:end-1)];
    p = -1; % counter to adjust for change in length
    for i = find(trBreak<minbreak(trialstage))
        p = p+1;
        A = trialdef{trialstage}(1,:);
        B = trialdef{trialstage}(2,:);
        % delete bad ones
        A(i-p) = [];
        B(i-p-1) = [];
        nmx(i-p) = [];
        trialdef{trialstage} = [A;B];
    end
    %%
    % Remove very small blocks
    trialdef{trialstage}(:,nmx<blockSize(trialstage)) = []; nmx(nmx<blockSize(trialstage)) = [];
end

% Bugfix for blocks in which there is leading null (code 1) before initial posture (code 3) (delete this
% trialdef row)
if  trialdef{1}(1,1)<trialdef{3}(1,1)
    trialdef{1}(:,1) = [];
end

% Switch between cases of subdivision
switch splitType
    
    case 'blockSplit'
        naming = {'Rest','Posture','Task'};
        restdef = trialdef{2}(:,1);
        postdef = trialdef{3}(:,1);
        taskdef = [trialdef{2}(2,1) trialdef{1}(1,1);
            trialdef{1}(2,1) trialdef{1}(1,2);
            trialdef{1}(2,2) trialdef{1}(1,3);
            trialdef{1}(2,3) trialdef{4}(end);
            ];
        
        trialdef = [{restdef} {postdef} {taskdef}];
        SMRtrialdef = [naming; trialdef]';
        SMRtrialdef = [SMRtrialdef repmat(trialdetails,size(trialdef,2),1)];
    case 'taskSplit'
        naming = {'MotorRest','MotorPrep','MotorExec','Hold'};
        motorRest = trialdef{1};
        motorPrep = trialdef{2};
        motorExec = trialdef{3};
        motorHold = trialdef{4};
        
        %         % Correct motorprep block one (extends into pretask region)
        %         motorPrep(1,1) = motorPrep(2,1)- fix(mean(diff(motorPrep(:,2:end))));
        
        trialdef = [{motorRest} {motorPrep} {motorExec} {motorHold}];
        SMRtrialdef = [naming; trialdef]';
        SMRtrialdef = [SMRtrialdef repmat(trialdetails,size(trialdef,2),1)];
        
    case 'condSplit'
        naming = {'1','2','3','4'};
        catCond = [trialdef{:}];

        % Script below runs through each condition marker
        for cond = 1:numel(naming)
            TD = trialdef{cond}(1,:);
            tmpDef = [];
            for p = 1:numel(TD)
                % find closest start of next trial
                trDist = TD(p)-catCond(1,:); % distance
                trDist(trDist>=0) = nan;
                if any(trDist)
                    [~,k] = max(trDist); % find minimum
                    tmpDef = [tmpDef [TD(p); catCond(1,k)]]; % construct trial
                else
                    tmpDef = [tmpDef [TD(p); size(ampCoder,2)]]; % construct trial
                end
            end
            condDef{cond} = tmpDef;
        end
        
        %% Construct the def
        trialdef = condDef;
        SMRtrialdef = [naming; trialdef]';
        SMRtrialdef = [SMRtrialdef repmat(trialdetails,size(trialdef,2),1)];
       
end
