function bridgeInds = bridgeCellSplits(Xinds,minwid)
if ~isempty(Xinds)
    p = 1;
    bridgeInds = {}; sp = [];
    bridgeInds{1} = Xinds{1};
    for i = 1:(numel(Xinds)-1)
        splitN = diff([Xinds{i}(end) Xinds{i+1}(1)]);
        if splitN< minwid
            bridgeInds{p} = [bridgeInds{p} Xinds{i+1}];
        else
            p = p+1; % or else split into new cell
            bridgeInds{p} = Xinds{i+1};
        end
        sp(i) = splitN;
    end
else
    bridgeInds = Xinds;
end
