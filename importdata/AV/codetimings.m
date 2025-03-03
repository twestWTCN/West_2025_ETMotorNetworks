function [decoder,SMRtrialdef] = codetimings(CVM,ampCoder,trialdetails,splitType,minblockSize,maxbreak,UNConds,maxblocksize)
if nargin<8
    maxblocksize = repmat(inf,1,size(CVM,1));
end
% ---- ampCoder = Coder coming from the amplifier of a movement/calibration trial (port #14)
% ---- CVM      = table extracted from trigger calibration measurement

% Try to decode using simple absolute differences
ampCoder(isnan(ampCoder)) = 0;
ampDeCoder = abs(ampCoder-CVM); % absolute Difference
[dError,decoder] = min(ampDeCoder,[],1);
dEx = 0.1.*min(diff(CVM));
decoder(dError>dEx) = nan;
% Now make a trialdef variable that you can use to segment the data in
% fieldtrip
condBank = cell(1,numel(CVM));
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
    
    
    %% Remove small breaks (minor collisions)- this basically eliminates trial definitions and concatanates them together
    % the break has to be less than maxbreak definition (i.e., max period of
    % time that you expect to occur in epoch.
     [trialdef{trialstage},nmx] = fixBreaks(trialdef{trialstage},maxbreak(trialstage),nmx);
    
    %%
    % Remove very small blocks
    trialdef{trialstage}(:,nmx<minblockSize(trialstage)) = []; nmx(nmx<minblockSize(trialstage)) = [];
    % Remove large blocks
    trialdef{trialstage}(:,nmx>maxblocksize(trialstage)) = []; nmx(nmx>maxblocksize(trialstage)) = [];
    
    condBank{trialstage} = [condBank{trialstage}  [trialdef{trialstage}(1,:);repmat(trialstage,1,size(trialdef{trialstage},2))]];
    
    %% If you want to check (can place above also
    % %     plot(ampCoder)
    % %     hold on
    % %     for i = 1:size(trialdef{trialstage},2)
    % %         plot(trialdef{trialstage}(:,i)',ampCoder(trialdef{trialstage}(:,i)'),'r','LineWidth',1.5)
    % %     end
end



% Switch between cases of subdivision
switch splitType
    
    case 'blockSplit'
        naming = {'Rest','Posture','Task'};
        restdef = trialdef{2}(:,1);
        postdef = trialdef{3}(:,1);
        taskdef = [trialdef{2}(2,1) trialdef{1}(1,2);
            trialdef{1}(2,2) trialdef{1}(1,3);
            trialdef{1}(2,3) trialdef{1}(1,4);
            trialdef{1}(2,4) trialdef{4}(end);
            ];
        
        trialdef = [{restdef} {postdef} {taskdef}];
        SMRtrialdef = [naming; trialdef]';
        SMRtrialdef = [SMRtrialdef repmat(trialdetails,size(trialdef,2),1)];
    case 'taskSplit'
        naming = {'MotorRest','MotorPost','MotorPrep','MotorExec','Hold'};
        motorRest = trialdef{1};
        motorPost = trialdef{2};
        motorPrep = trialdef{3};
        motorExec = trialdef{4};
        motorHold = trialdef{5};
        %         % Correct motorprep block one (extends into pretask region)
        %         motorPrep(1,1) = motorPrep(2,1)- fix(mean(diff(motorPrep(:,2:end))));
        
        trialdef = [{motorRest} {motorPost} {motorPrep} {motorExec} {motorHold}];
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
        
        X = horzcat(condBank{:});
        [XS,I] = sort(X(1,:));
        XMap = [XS; X(2,I)]; % this is derived from the decoder
        
        if numel(vertcat(UNConds{:,1})) == numel(X(2,I))
            if ~any(vertcat(UNConds{:,1})-X(2,I)')
                disp('Epoching matches Unity data! Adding behavioural info...')
                for cond = 1:numel(naming)
                    [C,IA] = intersect(XMap(1,:),SMRtrialdef{cond,2}(1,:));
                    SMRtrialdef{cond,5} = UNConds(IA,:)';
                end
            end
        else
            warning('Epoching doesnt match Unity data, try to rematch epochs...')
            % The script below finds the source of the mismatch
            XC = X(2,I); %[X(2,I) nan(1,sizediff)]';
            UNC = vertcat(UNConds{:,1});
            matchInd = 1; IC = [];
            for trl = 1:numel(XC)
                
                while XC(trl) ~= UNC(matchInd)
                    matchInd = matchInd+1;
                end
                IC(trl) = matchInd;
                matchInd = matchInd+1;
            end
            % and now do as above
            for cond = 1:numel(naming)
                [C,IA] = intersect(XMap(1,:),SMRtrialdef{cond,2}(1,:));
                SMRtrialdef{cond,5} = UNConds(IC(IA),:)';
            end
            disp('Successfully paired up the epochs with the Unity Data!')
        end
end
end
