function [pre,post] = getPrePostEpochDetails(part)  
% times of pre and post definitons
if part == 2 % post onset
    pre =  [-1.3 -0.500];
    post = [0 1];
elseif part == 3 % cue onset
    pre =  [-1 0];
    post = [0.1 1.1]; % to account for eyeblink
elseif part == 4 % move onset
    pre =  [-1.3 -0.300];
    post = [0.3 1.3];
elseif part == 5 %hold onset
    pre =  [-1 0];
    post = [1 2];
end