function [XData,YData,ZData] = setFTLengthToMin(XData,YData,ZData)
if nargin <3
    L = [size(XData.trial{1},2) size(YData.trial{1},2)];
else
    L = [size(XData.trial{1},2) size(YData.trial{1},2) size(ZData.trial{1},2)];
end

fsamp = XData.fsample;
L = min(L); % use the smallest
tvec = linspace(0,L/fsamp,L);

% Truncate
XData.trial{1} = XData.trial{1}(:,1:L);
YData.trial{1} = YData.trial{1}(:,1:L);
if nargin ==3
    ZData.trial{1} = ZData.trial{1}(:,1:L);
end
% Set Time
XData.time{1} = tvec;
YData.time{1} = tvec;
if nargin ==3
    ZData.time{1} = tvec;
end
% SampleInf
XData.sampleinfo = [0 L];
YData.sampleinfo = [0 L];
if nargin ==3
    ZData.sampleinfo = [0 L];
end