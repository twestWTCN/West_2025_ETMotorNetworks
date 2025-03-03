function [trialdef,nmx] = fixBreaks(trialdef,maxbreak,nmx)
    %% Remove small breaks (minor collisions)- this basically eliminates trial definitions and concatanates them together
    % the break has to be less than maxbreak definition (i.e., max period of
    % time that you expect to occur in epoch.
    trBreak = [inf trialdef(2,2:end)-trialdef(2,1:end-1)];
    p = -1; % counter to adjust for change in length
    for i = find(trBreak<maxbreak)
        p = p+1;
        A = trialdef(1,:);
        B = trialdef(2,:);
        % delete bad ones
        A(i-p) = [];
        B(i-p-1) = [];
        nmx(i-p) = [];
        trialdef = [A;B];
    end