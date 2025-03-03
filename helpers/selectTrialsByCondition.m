function SMRtrialdef = selectTrialsByCondition(SMRtrialdef,condtrialdef,condDeets)
% This function selects trials of "SMRtrialdef" if they are included within
% the condition specific "condtrialdef"

for block = 1:size(SMRtrialdef)
    % This is getting the full indices of data covered by this particular
    % block (e.g. condition)
    clist = [];
    for C = 1:size(condtrialdef,2)
        clist = [clist condtrialdef(1,C):condtrialdef(2,C)];
    end % this makes the complete list
    
    % This now checks for overlap between each reach and the main block
    ovlp = [];
    for trl = 1:size(SMRtrialdef{block,2},2)
        L =  SMRtrialdef{block,2}(1,trl): SMRtrialdef{block,2}(2,trl);
        ovlp(trl) = numel(intersect(L,clist))/numel(L); % find the overlap between
    end
    SMRtrialdef{block,2} = SMRtrialdef{block,2}(:,ovlp>0.9); % select reachs that overlap 90% with the condition marker
    
    % Now run through again and find which reach matches which Unity set
    condID = [];
    for trl = 1:size(SMRtrialdef{block,2},2)
        L = SMRtrialdef{block,2}(1,trl):SMRtrialdef{block,2}(2,trl);
        ovlp = [];
        for condtrl = 1:size(condtrialdef,2)
            C =  condtrialdef(1,condtrl):condtrialdef(2,condtrl);
            ovlp(condtrl) = numel(intersect(L,C))/numel(L);
        end
        [~,condID(trl)] =  max(ovlp);
    end
             
    SMRtrialdef{block,5} = condDeets(:,condID);
end